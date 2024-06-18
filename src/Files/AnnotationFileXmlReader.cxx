
/*LICENSE_START*/
/*
 *  Copyright (C) 2015 Washington University School of Medicine
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*LICENSE_END*/

#include <QFile>
#include <QXmlStreamAttributes>
#include <QXmlStreamReader>

#define __ANNOTATION_FILE_XML_READER_DECLARE__
#include "AnnotationFileXmlReader.h"
#undef __ANNOTATION_FILE_XML_READER_DECLARE__

#include "AnnotationBox.h"
#include "AnnotationCoordinate.h"
#include "AnnotationFile.h"
#include "AnnotationGroup.h"
#include "AnnotationImage.h"
#include "AnnotationLine.h"
#include "AnnotationMetaData.h"
#include "AnnotationOval.h"
#include "AnnotationPolygon.h"
#include "AnnotationPolyLine.h"
#include "AnnotationPolyhedron.h"
#include "AnnotationPercentSizeText.h"
#include "AnnotationPointSizeText.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "DataFileException.h"
#include "GiftiXmlElements.h"
#include "XmlStreamReaderHelper.h"

using namespace caret;


    
/**
 * \class caret::AnnotationFileXmlReader 
 * \brief Read Annotation File from XML format.
 * \ingroup Files
 */

/**
 * Constructor.
 */
AnnotationFileXmlReader::AnnotationFileXmlReader()
: AnnotationFileXmlFormatBase()
{
    m_stream.grabNew(new QXmlStreamReader());
    m_streamHelper.grabNew(NULL);
}

/**
 * Destructor.
 */
AnnotationFileXmlReader::~AnnotationFileXmlReader()
{
}

/**
 * Read the given annotation file in XML format.
 *
 * @param filename
 *    Name of file.
 * @param annotationFile
 *    Read into this annotation file.
 * @throws
 *    DataFileException if there is an error reading the file.
 */
void
AnnotationFileXmlReader::readFile(const QString& filename,
                                  AnnotationFile* annotationFile)
{
    m_fileVersionNumber = -1;
    
    CaretAssert(annotationFile);
    m_filename = filename;
    setAnnotationFileDirectory(filename);

    /*
     * Open the file
     */
    QFile file(m_filename);
    if ( ! file.open(QFile::ReadOnly)) {
        const QString msg("Error opening for reading."
                          + file.errorString());
        throw DataFileException(m_filename,
                                msg);
    }
    
    /*
     * Create an XML stream writer
     */
    m_stream.grabNew(new QXmlStreamReader(&file));
    
    readFileContentFromXmlStreamReader(filename,
                                       annotationFile);
    
    file.close();
    
    if (m_stream->hasError()) {
        m_streamHelper->throwDataFileException("There was an error reading the annotation file in XML format (reported by QXmlStreamReader): "
                                + m_stream->errorString());
    }
}

/**
 * Read the given annotation file contained in a string in XML format.
 *
 * @param fileInString
 *    String containing the file.
 * @param fileNameForRelativePaths
 *    Filename used when relative paths are stored in the annotation file
 * @param annotationFile
 *    Read into this annotation file.
 * @throws
 *    DataFileException if there is an error reading the file.
 */
void
AnnotationFileXmlReader::readFileFromString(const QString& fileInString,
                                            const AString& fileNameForRelativePaths,
                                            AnnotationFile* annotationFile)
{
    /*
     * Create an XML stream writer
     */
    m_stream.grabNew(new QXmlStreamReader(fileInString));
    setAnnotationFileDirectory(fileNameForRelativePaths);
    
    if (annotationFile != NULL) {
        m_filename = annotationFile->getFileName();
    }
    readFileContentFromXmlStreamReader("AnnotationsInSceneFile",
                                       annotationFile);
    
    if (m_stream->hasError()) {
        m_streamHelper->throwDataFileException("There was an error reading the annotation file in XML format (reported by QXmlStreamReader): "
                                               + m_stream->errorString());
    }
}

/**
 * Read content of the file from the XML stream reader.
 *
 * @param filename
 *    Name of the file.
 * @param annotationFile
 *    Read into this annotation file.
 * @throws
 *    DataFileException if there is an error reading the file.
 */
void
AnnotationFileXmlReader::readFileContentFromXmlStreamReader(const QString& filename,
                                        AnnotationFile* annotationFile)
{
    /*
     * Create the helper for reading XML
     */
    if (m_streamHelper != NULL) {
        m_streamHelper.grabNew(NULL);
    }
    m_streamHelper.grabNew(new XmlStreamReaderHelper(filename,
                                                     m_stream));
    
    if (m_stream->atEnd()) {
        m_streamHelper->throwDataFileException("Error reading.  File appears to have no XML content.");
    }
    
    const bool fileElementValid = m_stream->readNextStartElement();
    if ( ! fileElementValid) {
        m_streamHelper->throwDataFileException("Appears to have no XML elements.");
    }
    
    const auto fileElementName = m_stream->name();
    if (fileElementName != ELEMENT_ANNOTATION_FILE) {
        m_streamHelper->throwDataFileException("First element is "
                                               + fileElementName.toString()
                                               + " but should be "
                                               + ELEMENT_ANNOTATION_FILE);
    }
    
    QXmlStreamAttributes fileAttributes = m_stream->attributes();
    const auto versionText = fileAttributes.value(ATTRIBUTE_VERSION);
    if (versionText.isEmpty()) {
        m_streamHelper->throwDataFileException("Version attribute ("
                                               + ATTRIBUTE_VERSION
                                               + ") is missing from the file element "
                                               + ELEMENT_ANNOTATION_FILE);
    }
    
    m_fileVersionNumber  = versionText.toString().toInt();
    if (m_fileVersionNumber == XML_VERSION_ONE) {
        readVersionOne(annotationFile);
    }
    else if (m_fileVersionNumber == XML_VERSION_TWO) {
        readVersionTwo(annotationFile);
    }
    else if (m_fileVersionNumber == XML_VERSION_THREE) {
        /* 
         * NOTE: version 3 added new coordinate space "SPACER "
         * and otherwise is the same as version 2 format
         */
        readVersionTwo(annotationFile);
    }
    else {
        m_streamHelper->throwDataFileException("File version number "
                                               + versionText.toString()
                                               + " is not supported by this version of the software.");
    }
}

/**
 * Read a version one Annotation XML file.
 *
 * @param annotationFile
 *     Add annotations to this file.
 */
void
AnnotationFileXmlReader::readVersionOne(AnnotationFile* annotationFile)
{
    while (m_stream->readNextStartElement()) {
        bool skipCurrentElementFlag = true;
        
        const QString elementName = m_stream->name().toString();
        if (elementName == GiftiXmlElements::TAG_METADATA) {
            m_streamHelper->readMetaData(annotationFile->getFileMetaData());
            skipCurrentElementFlag = false;
        }
        else if (elementName == ELEMENT_BOX) {
            CaretPointer<AnnotationBox> annotation(new AnnotationBox(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_BOX,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_IMAGE) {
            CaretPointer<AnnotationImage> annotation(new AnnotationImage(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_IMAGE,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_LINE) {
            CaretPointer<AnnotationLine> annotation(new AnnotationLine(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readTwoCoordinateAnnotation(annotationFile,
                                        ELEMENT_LINE,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_OVAL) {
            CaretPointer<AnnotationOval> annotation(new AnnotationOval(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_OVAL,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_PERCENT_SIZE_TEXT) {
            CaretPointer<AnnotationText> annotation(new AnnotationPercentSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_PERCENT_SIZE_TEXT,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_PERCENT_WIDTH_SIZE_TEXT) {
            CaretPointer<AnnotationText> annotation(new AnnotationPercentSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL,
                                                                                  AnnotationTextFontSizeTypeEnum::PERCENTAGE_OF_VIEWPORT_WIDTH));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_PERCENT_WIDTH_SIZE_TEXT,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_POINT_SIZE_TEXT) {
            CaretPointer<AnnotationText> annotation(new AnnotationPointSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_POINT_SIZE_TEXT,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_TEXT_OBSOLETE) {
            CaretPointer<AnnotationText> annotation(new AnnotationPercentSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_TEXT_OBSOLETE,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else {
            m_streamHelper->throwDataFileException("Unexpected XML element "
                                   + elementName);
        }

        /*
         * These elements have no other child elements so move on
         */
        if (skipCurrentElementFlag) {
            m_stream->skipCurrentElement();
        }
    }
}

/**
 * Read a version two Annotation XML file.
 *
 * @param annotationFile
 *     Add annotations to this file.
 */
void
AnnotationFileXmlReader::readVersionTwo(AnnotationFile* annotationFile)
{
    while (m_stream->readNextStartElement()) {
        bool skipCurrentElementFlag = true;
        
        const QString elementName = m_stream->name().toString();
        
        if (elementName == GiftiXmlElements::TAG_METADATA) {
            m_streamHelper->readMetaData(annotationFile->getFileMetaData());
            skipCurrentElementFlag = false;
        }
        else if (elementName == ELEMENT_GROUP) {
            readGroup(annotationFile);
            skipCurrentElementFlag = false;
        }
        else {
            m_streamHelper->throwDataFileException("Unexpected XML element "
                                                   + elementName);
        }
        
        /*
         * These elements have no other child elements so move on
         */
        if (skipCurrentElementFlag) {
            m_stream->skipCurrentElement();
        }
    }
}

/**
 * Read the next start element which should be a coordinate
 * with the given element name.
 *
 * @param coordinateElementName
 *     Element name to read for coordinate
 * @param coordinate
 *     Coordinate whose data is set with data read
 * @param coordinateSpace
 *     Coordinate space of annotation
 *  @param readStartElementFlag
 *     If true, read the start element and skipt to next element when done
 * @throw
 *     DataFileException
 */
void
AnnotationFileXmlReader::readCoordinate(const QString& coordinateElementName,
                                        AnnotationCoordinate* coordinate,
                                        const AnnotationCoordinateSpaceEnum::Enum coordinateSpace,
                                        const bool readStartElementFlag)
{
    CaretAssert(coordinate);
    
    if (readStartElementFlag) {
        const bool elementValid = m_stream->readNextStartElement();
        if ( ! elementValid) {
            m_streamHelper->throwDataFileException("Failed to read element "
                                                   + coordinateElementName);
        }
        
        if (m_stream->name() != coordinateElementName) {
            m_streamHelper->throwDataFileException("Expected elment "
                                                   + coordinateElementName
                                                   + " but read element "
                                                   + m_stream->name().toString());
        }
    }
    
    const QXmlStreamAttributes attributes = m_stream->attributes();

    /*
     * XYZ coordinate
     */
    const float xyz[3] = {
        m_streamHelper->getRequiredAttributeFloatValue(attributes, coordinateElementName, ATTRIBUTE_COORD_X),
        m_streamHelper->getRequiredAttributeFloatValue(attributes, coordinateElementName, ATTRIBUTE_COORD_Y),
        m_streamHelper->getRequiredAttributeFloatValue(attributes, coordinateElementName, ATTRIBUTE_COORD_Z)
    };
    coordinate->setXYZ(xyz);

    /*
     * Surface coordinate
     */
    const int32_t numberOfNodes = m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                 coordinateElementName,
                                                                 ATTRIBUTE_COORD_SURFACE_NUMBER_OF_NODES);
    const int32_t nodeIndex = m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                 coordinateElementName,
                                                                 ATTRIBUTE_COORD_SURFACE_NODE_INDEX);
    const QString structureValueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                coordinateElementName,
                                                                ATTRIBUTE_COORD_SURFACE_STRUCTURE);
    const float offsetDistance = m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                                                coordinateElementName,
                                                                                ATTRIBUTE_COORD_SURFACE_NODE_OFFSET);

    const QString offsetVectorString = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                       coordinateElementName,
                                                                                       ATTRIBUTE_COORD_SURFACE_NODE_OFFSET_VECTOR_TYPE,
                                                                                       "CENTROID_THRU_VERTEX");
    bool structureValid = false;
    StructureEnum::Enum structure = StructureEnum::fromName(structureValueString,
                                                          &structureValid);

    bool offsetVectorValid = false;
    AnnotationSurfaceOffsetVectorTypeEnum::Enum offsetVector = AnnotationSurfaceOffsetVectorTypeEnum::fromName(offsetVectorString,
                                                                                                               &offsetVectorValid);


    if ( ! offsetVectorValid) {
        offsetVector = AnnotationSurfaceOffsetVectorTypeEnum::CENTROID_THRU_VERTEX;
    }

    if (structureValid) {
        coordinate->setSurfaceSpace(structure,
                                    numberOfNodes,
                                    nodeIndex,
                                    offsetDistance,
                                    offsetVector);
    }
    else {
        m_streamHelper->throwDataFileException("Invalid value "
                               + structureValueString
                               + " for attribute "
                               + ATTRIBUTE_COORD_SURFACE_STRUCTURE);
    }

    switch (coordinateSpace) {
        case AnnotationCoordinateSpaceEnum::MEDIA_FILE_NAME_AND_PIXEL:
        {
            if (m_stream->readNextStartElement()) {
                if (m_stream->name() == ELEMENT_COORDINATE_MEDIA_FILE_NAME) {
                    const QString mediaFileName = m_stream->readElementText(QXmlStreamReader::ErrorOnUnexpectedElement);
                    coordinate->setMediaFileName(mediaFileName);
                }
                else {
                    m_streamHelper->throwDataFileException("Expected elment "
                                                           + ELEMENT_COORDINATE_MEDIA_FILE_NAME
                                                           + " but read element "
                                                           + m_stream->name().toString());
                }
            }
            else {
                m_streamHelper->throwDataFileException("Failed to coordinate child element "
                                                       + ELEMENT_COORDINATE_MEDIA_FILE_NAME);
            }
        }
            break;
        case AnnotationCoordinateSpaceEnum::CHART:
            break;
        case AnnotationCoordinateSpaceEnum::HISTOLOGY:
        {
            if (m_stream->readNextStartElement()) {
                if (m_stream->name() == ELEMENT_COORDINATE_HISTOLOGY_SPACE_KEY) {
                    const AString encodedHistologyKey(m_stream->readElementText(QXmlStreamReader::ErrorOnUnexpectedElement));
                    HistologySpaceKey hsk;
                    if (hsk.setFromEncodedString(getAnnotationFileDirectory(),
                                                 encodedHistologyKey)) {
                        coordinate->setHistologySpaceKey(hsk);
                    }
                    else {
                        m_streamHelper->throwDataFileException(("Failed to decode \""
                                                                + encodedHistologyKey
                                                                + " into a HistologyKey"));
                    }
                }
                else {
                    m_streamHelper->throwDataFileException("Expected elment "
                                                           + ELEMENT_COORDINATE_HISTOLOGY_SPACE_KEY
                                                           + " but read element "
                                                           + m_stream->name().toString());
                }
            }
            else {
                m_streamHelper->throwDataFileException("Failed to find coordinate child element "
                                                       + ELEMENT_COORDINATE_HISTOLOGY_SPACE_KEY);
            }
        }
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
        case AnnotationCoordinateSpaceEnum::SURFACE:
        case AnnotationCoordinateSpaceEnum::TAB:
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
        case AnnotationCoordinateSpaceEnum::WINDOW:
            break;
    }

    if (readStartElementFlag) {
        /*
         * No other child elements so move on
         */
        m_stream->skipCurrentElement();
    }
}

/**
 * Read the attributes common to all annotation elements.
 *
 * @param annotation
 *     The annotation.
 * @param annotationElementName
 *     Name of element of annotation.
 * @param attributes
 *     The XML attributes.
 * @throw
 *     DataFileException if there is an error reading the attributes.
 */
void
AnnotationFileXmlReader::readAnnotationAttributes(Annotation* annotation,
                                                  const QString& annotationElementName,
                                                  const QXmlStreamAttributes& attributes)
{
    CaretAssert(annotation);
    
    {
        /*
         * Coordinate space
         */
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    annotationElementName,
                                                                    ATTRIBUTE_COORDINATE_SPACE);
        bool valid = false;
        AnnotationCoordinateSpaceEnum::Enum value = AnnotationCoordinateSpaceEnum::fromName(valueString,
                                                                                            &valid);
        if (valid) {
            annotation->setCoordinateSpace(value);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_COORDINATE_SPACE);
        }
    }
    
    {
        /*
         * Background color
         */
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    annotationElementName,
                                                                    ATTRIBUTE_BACKGROUND_CARET_COLOR);
        bool valid = false;
        CaretColorEnum::Enum value = CaretColorEnum::fromName(valueString,
                                                              &valid);
        if (valid) {
            annotation->setBackgroundColor(value);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_BACKGROUND_CARET_COLOR);
        }
    }
    
    {
        /*
         * Background custom color
         */
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    annotationElementName,
                                                                    ATTRIBUTE_BACKGROUND_CUSTOM_RGBA);
        std::vector<float> rgba;
        AString::toNumbers(valueString, rgba);
        if (rgba.size() == 4) {
            annotation->setCustomBackgroundColor(&rgba[0]);
        }
        else {
            m_streamHelper->throwDataFileException(ATTRIBUTE_BACKGROUND_CUSTOM_RGBA
                                   + " must contain 4 elements but "
                                   + valueString
                                   + " contains "
                                   + QString::number(rgba.size())
                                   + " elements");
        }
    }
    
    {
        /*
         * Foreground color
         */
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    annotationElementName,
                                                                    ATTRIBUTE_FOREGROUND_CARET_COLOR);
        bool valid = false;
        CaretColorEnum::Enum value = CaretColorEnum::fromName(valueString,
                                                              &valid);
        if (valid) {
            annotation->setLineColor(value);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_FOREGROUND_CARET_COLOR);
        }
    }
    
    {
        /*
         * Foreground custom color
         */
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    annotationElementName,
                                                                    ATTRIBUTE_FOREGROUND_CUSTOM_RGBA);
        std::vector<float> rgba;
        AString::toNumbers(valueString, rgba);
        if (rgba.size() == 4) {
            annotation->setCustomLineColor(&rgba[0]);
        }
        else {
            m_streamHelper->throwDataFileException(ATTRIBUTE_FOREGROUND_CUSTOM_RGBA
                                   + " must contain 4 elements but "
                                   + valueString
                                   + " contains "
                                   + QString::number(rgba.size())
                                   + " elements");
        }
    }
    
    /*
     * Foreground line width
     */
    annotation->setLineWidthPixelsObsolete(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                                      annotationElementName,
                                                                      ATTRIBUTE_FOREGROUND_LINE_WIDTH_PIXELS));

    /*
     * Line width percentage added on July 26, 2017.
     * A "negative value" indicates it is missing and the percentage width
     * may be set by graphics code that attempts to set the percentage using
     * the obsolete pixel width and viewport height.
     */
    annotation->setLineWidthPercentage(m_streamHelper->getOptionalAttributeFloatValue(attributes,
                                                                                    annotationElementName,
                                                                                    ATTRIBUTE_FOREGROUND_LINE_WIDTH_PERCENTAGE,
                                                                                    -1.0f));
    /*
     * Tab Index
     */
    annotation->setTabIndex(m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                         annotationElementName,
                                                         ATTRIBUTE_TAB_INDEX));
    /*
     * Window Index
     */
    annotation->setWindowIndex(m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                            annotationElementName,
                                                            ATTRIBUTE_WINDOW_INDEX));
    /*
     * Spacer Tab Index added as part of WB-668
     */
    const AString spacerTabText = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                 annotationElementName,
                                                                                 ATTRIBUTE_SPACER_TAB_INDEX,
                                                                                 "");
    SpacerTabIndex spacerTabIndex;
    if ( ! spacerTabText.isEmpty()) {
        spacerTabIndex.setFromXmlAttributeText(spacerTabText);
    }
    annotation->setSpacerTabIndex(spacerTabIndex);
    
    /*
     * Unique Key
     */
    if (m_fileVersionNumber >= XML_VERSION_TWO) {
        annotation->setUniqueKey(m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                              annotationElementName,
                                                                              ATTRIBUTE_UNIQUE_KEY));
    }
}

/**
 * Read a two coordinate annotation.
 *
 * @param annotationElementName
 *     Name of one-dimensional attribute.
 * @param annotation
 *     One-dimensional annotation that has its data read.
 */
void
AnnotationFileXmlReader::readTwoCoordinateAnnotation(AnnotationFile* annotationFile,
                                                     const QString& annotationElementName,
                                                     AnnotationTwoCoordinateShape* annotation)
{
    CaretAssert(annotation);

    const QXmlStreamAttributes attributes = m_stream->attributes();
    
    readAnnotationAttributes(annotation,
                             annotationElementName,
                             attributes);
    
    AnnotationLine* line = dynamic_cast<AnnotationLine*>(annotation);
    if (line != NULL) {
        line->setDisplayEndArrow(m_streamHelper->getOptionalAttributeBoolValue(attributes,
                                                                               annotationElementName,
                                                                               ATTRIBUTE_LINE_END_ARROW,
                                                                               false));
        line->setDisplayStartArrow(m_streamHelper->getOptionalAttributeBoolValue(attributes,
                                                                                 annotationElementName,
                                                                                 ATTRIBUTE_LINE_START_ARROW,
                                                                                 false));
    }
    
    bool done(false);
    while ( ! done) {
        const QXmlStreamReader::TokenType tokenType(m_stream->readNext());
        if (m_stream->atEnd()) {
            done = true;
        }
        else {
            const QString elementName(m_stream->name().toString());
            switch (tokenType) {
                case QXmlStreamReader::StartElement:
                    if (elementName == ELEMENT_COORDINATE_ONE) {
                        /*
                         * Read the coordinate
                         */
                        const bool readStartElementFlag(false);
                        readCoordinate(ELEMENT_COORDINATE_ONE,
                                       annotation->getStartCoordinate(),
                                       annotation->getCoordinateSpace(),
                                       readStartElementFlag);
                        m_stream->skipCurrentElement();
                    }
                    else if (elementName == ELEMENT_COORDINATE_TWO) {
                        /*
                         * Read the coordinate
                         */
                        const bool readStartElementFlag(false);
                        readCoordinate(ELEMENT_COORDINATE_TWO,
                                       annotation->getEndCoordinate(),
                                       annotation->getCoordinateSpace(),
                                       readStartElementFlag);
                        m_stream->skipCurrentElement();
                    }
                    else if (elementName == GiftiXmlElements::TAG_METADATA) {
                        /*
                         * Reads through the metadata end element so
                         * DO NOT skip to next element
                         */
                        m_streamHelper->readMetaData(annotation->getMetaData());
                    }
                    else {
                        /*
                         * Issue warning (instead of fatal error) if unrecognized element found.
                         * Will skip over remainder of element later in code.
                         */
                        annotationFile->addFileReadWarning("Unrecognized element \""
                                                           + elementName
                                                           + "\" and content ignored.  Updating Workbench "
                                                           "may correct this problem.");
                        m_stream->skipCurrentElement();
                    }
                    break;
                case QXmlStreamReader::EndElement:
                    if (elementName == annotationElementName) {
                        done = true;
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

/**
 * Read a multi coordinate annotation.
 *
 * @param annotationFile,
       Annotation file being read
 * @param annotationElementName
 *     Name of one-dimensional attribute.
 * @param annotation
 *     One-dimensional annotation that has its data read.
 */
void
AnnotationFileXmlReader::readMultiCoordinateAnnotation(AnnotationFile* annotationFile,
                                                       const QString& annotationElementName,
                                                       AnnotationMultiCoordinateShape* annotation)
{
    CaretAssert(annotation);
    
    const QXmlStreamAttributes attributes = m_stream->attributes();
    
    readAnnotationAttributes(annotation,
                             annotationElementName,
                             attributes);
    
    bool done(false);
    while ( ! done) {
        const QXmlStreamReader::TokenType tokenType(m_stream->readNext());
        if (m_stream->atEnd()) {
            done = true;
        }
        else {
            const QString elementName(m_stream->name().toString());
            switch (tokenType) {
                case QXmlStreamReader::StartElement:
                    if (elementName == ELEMENT_COORDINATE_LIST) {
                        const QXmlStreamAttributes coordListAtts(m_stream->attributes());
                        const int32_t numberOfCoordiantes = m_streamHelper->getRequiredAttributeIntValue(coordListAtts,
                                                                                                         ELEMENT_COORDINATE_LIST,
                                                                                                         ATTRIBUTE_COORDINATE_LIST_COUNT);
                        for (int32_t i = 0; i < numberOfCoordiantes; i++) {
                            AnnotationCoordinate* ac = new AnnotationCoordinate(annotation->m_attributeDefaultType);
                            const bool readStartElementFlag(true);
                            readCoordinate(ELEMENT_COORDINATE,
                                           ac,
                                           annotation->getCoordinateSpace(),
                                           readStartElementFlag);
                            annotation->addCoordinate(ac);
                        }
                        
                        m_stream->skipCurrentElement();
                    }
                    else if (elementName == GiftiXmlElements::TAG_METADATA) {
                        /*
                         * Reads through the metadata end element so
                         * DO NOT skip to next element
                         */
                        m_streamHelper->readMetaData(annotation->getMetaData());
                    }
                    else {
                        /*
                         * Issue warning (instead of fatal error) if unrecognized element found.
                         * Will skip over remainder of element later in code.
                         */
                        annotationFile->addFileReadWarning("Unrecognized element \""
                                                           + elementName
                                                           + "\" and content ignored.  Updating Workbench "
                                                           "may correct this problem.");
                        m_stream->skipCurrentElement();
                    }
                    break;
                case QXmlStreamReader::EndElement:
                    if (elementName == annotationElementName) {
                        done = true;
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

/**
 * Read a multi paired coordinate annotation.
 *
 * @param annotationFile
 *     Annnotation file being read
 * @param annotationElementName
 *     Name of one-dimensional attribute.
 * @param annotation
 *     One-dimensional annotation that has its data read.
 */
void
AnnotationFileXmlReader::readMultiPairedCoordinateAnnotation(AnnotationFile* annotationFile,
                                                             const QString& annotationElementName,
                                                             AnnotationMultiPairedCoordinateShape* annotation)
{
    CaretAssert(annotation);
    
    const QXmlStreamAttributes attributes = m_stream->attributes();
    
    readAnnotationAttributes(annotation,
                             annotationElementName,
                             attributes);

    AnnotationPolyhedron* polyhedron(annotation->castToPolyhedron());

    bool done(false);
    while ( ! done) {
        const QXmlStreamReader::TokenType tokenType(m_stream->readNext());
        if (m_stream->atEnd()) {
            done = true;
        }
        else {
            const QString elementName(m_stream->name().toString());
            switch (tokenType) {
                case QXmlStreamReader::StartElement:
                    if (elementName == ELEMENT_COORDINATE_LIST) {
                        const QXmlStreamAttributes coordListAtts(m_stream->attributes());
                        const int32_t numberOfCoordiantes = m_streamHelper->getRequiredAttributeIntValue(coordListAtts,
                                                                                                         ELEMENT_COORDINATE_LIST,
                                                                                                         ATTRIBUTE_COORDINATE_LIST_COUNT);
                        for (int32_t i = 0; i < numberOfCoordiantes; i++) {
                            AnnotationCoordinate* ac = new AnnotationCoordinate(annotation->m_attributeDefaultType);
                            const bool readStartElementFlag(true);
                            readCoordinate(ELEMENT_COORDINATE,
                                           ac,
                                           annotation->getCoordinateSpace(),
                                           readStartElementFlag);
                            annotation->addCoordinate(ac);
                        }
                        
                        m_stream->skipCurrentElement();
                        
                        if (polyhedron != NULL) {
                            /*
                             * XYZ for names is not in older polyhedrons.  This call
                             * will initialize name XYZ to center of coords on each plane
                             */
                            polyhedron->resetPlaneOneTwoNameStereotaxicXYZ();
                        }
                    }
                    else if (elementName == GiftiXmlElements::TAG_METADATA) {
                        /*
                         * Reads through the metadata end element so
                         * DO NOT skip to next element
                         */
                        m_streamHelper->readMetaData(annotation->getMetaData());
                    }
                    else if (elementName == ELEMENT_POLYHEDRON_DATA) {
                        CaretAssert(polyhedron);
                        
                        const QXmlStreamAttributes polyAtts(m_stream->attributes());

                        const AString planeOneString(m_streamHelper->getOptionalAttributeStringValue(polyAtts,
                                                                                                     ELEMENT_POLYHEDRON_DATA,
                                                                                                     ATTRIBUTE_PLANE_ONE,
                                                                                                     ""));
                        Plane planeOne;
                        if ( ! planeOneString.isEmpty()) {
                            planeOne = Plane::fromFormattedString(planeOneString);
                        }

                        const AString planeTwoString(m_streamHelper->getOptionalAttributeStringValue(polyAtts,
                                                                                                     ELEMENT_POLYHEDRON_DATA,
                                                                                                     ATTRIBUTE_PLANE_TWO,
                                                                                                     ""));
                        Plane planeTwo;
                        if ( ! planeTwoString.isEmpty()) {
                            planeTwo = Plane::fromFormattedString(planeTwoString);
                        }
                        
                        polyhedron->setFromFileReading(planeOne,
                                                       planeTwo);
                        


                        const AString planeOneNameXyzString(m_streamHelper->getOptionalAttributeStringValue(polyAtts,
                                                                                                            ELEMENT_POLYHEDRON_DATA,
                                                                                                            ATTRIBUTE_PLANE_ONE_NAME_XYZ,
                                                                                                            ""));
                        if ( ! planeOneNameXyzString.isEmpty()) {
                            bool validFlag(false);
                            const Vector3D xyz(Vector3D::fromString(planeOneNameXyzString,
                                                                    &validFlag));
                            if (validFlag) {
                                polyhedron->setPlaneOneNameStereotaxicXYZ(xyz);
                            }
                            else {
                                annotationFile->addFileReadWarning("Failed to convert \""
                                                                   + planeOneNameXyzString
                                                                   + "\" to a Vector3D");
                            }
                        }

                        const AString planeTwoNameXyzString(m_streamHelper->getOptionalAttributeStringValue(polyAtts,
                                                                                                            ELEMENT_POLYHEDRON_DATA,
                                                                                                            ATTRIBUTE_PLANE_TWO_NAME_XYZ,
                                                                                                            ""));
                        if ( ! planeTwoNameXyzString.isEmpty()) {
                            bool validFlag(false);
                            const Vector3D xyz(Vector3D::fromString(planeTwoNameXyzString,
                                                                    &validFlag));
                            if (validFlag) {
                                polyhedron->setPlaneTwoNameStereotaxicXYZ(xyz);
                            }
                            else {
                                annotationFile->addFileReadWarning("Failed to convert \""
                                                                   + planeTwoNameXyzString
                                                                   + "\" to a Vector3D");
                            }
                        }
                        
                        m_stream->skipCurrentElement();
                    }
                    else if (elementName == ELEMENT_FONT_ATTRIBUTES) {
                        AnnotationFontAttributesInterface* fontAttributesInterface(dynamic_cast<AnnotationFontAttributesInterface*>(annotation));
                        if (fontAttributesInterface != NULL) {
                            readFontAttibutes(fontAttributesInterface,
                                              elementName,
                                              m_stream->attributes());
                        }
                    }
                    else {
                        /*
                         * Issue warning (instead of fatal error) if unrecognized element found.
                         * Will skip over remainder of element later in code.
                         */
                        annotationFile->addFileReadWarning("Unrecognized element \""
                                                           + elementName
                                                           + "\" and content ignored.  Updating Workbench "
                                                           "may correct this problem.");
                        m_stream->skipCurrentElement();
                    }
                    break;
                case QXmlStreamReader::EndElement:
                    if (elementName == annotationElementName) {
                        done = true;
                    }
                    break;
                default:
                    break;
            }
        }
    }
    
    if (polyhedron != NULL) {
        /* Setting name xyz sets modified status so clear it */
        polyhedron->clearModified();
    }
}

/**
 * Read an annotation group.
 *
 * @param annotationFile
 *     File that is being read.
 */
void
AnnotationFileXmlReader::readGroup(AnnotationFile* annotationFile)
{
    const QXmlStreamAttributes attributes = m_stream->attributes();

    /*
     * Coordinate space
     */
    AnnotationCoordinateSpaceEnum::Enum coordSpace = AnnotationCoordinateSpaceEnum::VIEWPORT;
    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                                    ELEMENT_GROUP,
                                                                                    ATTRIBUTE_COORDINATE_SPACE);
        bool valid = false;
        coordSpace = AnnotationCoordinateSpaceEnum::fromName(valueString,
                                                                                            &valid);
        if ( ! valid) {
            m_streamHelper->throwDataFileException("While reading annotation group, invalid value "
                                                   + valueString
                                                   + " for attribute "
                                                   + ATTRIBUTE_COORDINATE_SPACE);
        }
    }

    /*
     * Group tyoe
     */
    AnnotationGroupTypeEnum::Enum groupType = AnnotationGroupTypeEnum::INVALID;
    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                                    ELEMENT_GROUP,
                                                                                    ATTRIBUTE_GROUP_TYPE);
        bool valid = false;
        groupType = AnnotationGroupTypeEnum::fromName(valueString,
                                                             &valid);
        if ( ! valid) {
            m_streamHelper->throwDataFileException("While reading annotation group, invalid value "
                                                   + valueString
                                                   + " for attribute "
                                                   + ATTRIBUTE_GROUP_TYPE);
        }
    }
    
    const int32_t tabOrWindowIndex = m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                                  ELEMENT_GROUP,
                                                                                  ATTRIBUTE_TAB_OR_WINDOW_INDEX);
    
    /*
     * Spacer Tab Index added as part of WB-668
     */
    const AString spacerTabText = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                  ELEMENT_GROUP,
                                                                                  ATTRIBUTE_SPACER_TAB_INDEX,
                                                                                  "");
    SpacerTabIndex spacerTabIndex;
    if ( ! spacerTabText.isEmpty()) {
        spacerTabIndex.setFromXmlAttributeText(spacerTabText);
    }
    
    const int32_t uniqueKey = m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                                  ELEMENT_GROUP,
                                                                                  ATTRIBUTE_UNIQUE_KEY);
    
    
    std::vector<Annotation*> annotations;
    AString mediaFileName;
    HistologySpaceKey histologySpaceKey;
    
    while (m_stream->readNextStartElement()) {
        const QString elementName = m_stream->name().toString();
        
        if (elementName == ELEMENT_COORDINATE_MEDIA_FILE_NAME) {
            mediaFileName = m_stream->readElementText(QXmlStreamReader::ErrorOnUnexpectedElement);
        }
        else if (elementName == ELEMENT_COORDINATE_HISTOLOGY_SPACE_KEY) {
            const AString encodedHistologyKey(m_stream->readElementText(QXmlStreamReader::ErrorOnUnexpectedElement));
            if ( ! histologySpaceKey.setFromEncodedString(getAnnotationFileDirectory(),
                                                          encodedHistologyKey)) {
                m_streamHelper->throwDataFileException(("Failed to decode \""
                                                        + encodedHistologyKey
                                                        + " into a HistologyKey"));
            }
        }
        else if (elementName == ELEMENT_BOX) {
            CaretPointer<AnnotationBox> annotation(new AnnotationBox(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_BOX,
                                         annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_IMAGE) {
            CaretPointer<AnnotationImage> annotation(new AnnotationImage(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_IMAGE,
                                         annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_LINE) {
            CaretPointer<AnnotationLine> annotation(new AnnotationLine(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readTwoCoordinateAnnotation(annotationFile,
                                        ELEMENT_LINE,
                                         annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_OVAL) {
            CaretPointer<AnnotationOval> annotation(new AnnotationOval(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_OVAL,
                                         annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_POLYGON) {
            CaretPointer<AnnotationPolygon> annotation(new AnnotationPolygon(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readMultiCoordinateAnnotation(annotationFile,
                                          ELEMENT_POLYGON,
                                           annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_POLY_LINE) {
            CaretPointer<AnnotationPolyLine> annotation(new AnnotationPolyLine(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readMultiCoordinateAnnotation(annotationFile,
                                          ELEMENT_POLY_LINE,
                                          annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_POLYHEDRON) {
            CaretPointer<AnnotationPolyhedron> annotation(new AnnotationPolyhedron(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readMultiPairedCoordinateAnnotation(annotationFile,
                                                ELEMENT_POLYHEDRON,
                                                annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_PERCENT_SIZE_TEXT) {
            CaretPointer<AnnotationText> annotation(new AnnotationPercentSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_PERCENT_SIZE_TEXT,
                                         annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_PERCENT_WIDTH_SIZE_TEXT) {
            CaretPointer<AnnotationText> annotation(new AnnotationPercentSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_PERCENT_WIDTH_SIZE_TEXT,
                                         annotation);
            annotationFile->addAnnotationDuringFileVersionOneReading(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_POINT_SIZE_TEXT) {
            CaretPointer<AnnotationText> annotation(new AnnotationPointSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_POINT_SIZE_TEXT,
                                         annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else if (elementName == ELEMENT_TEXT_OBSOLETE) {
            CaretPointer<AnnotationText> annotation(new AnnotationPercentSizeText(AnnotationAttributesDefaultTypeEnum::NORMAL));
            readOneCoordinateAnnotation(annotationFile,
                                        ELEMENT_TEXT_OBSOLETE,
                                        annotation);
            annotations.push_back(annotation.releasePointer());
        }
        else {
            /*
             * Issue warning (instead of fatal error) if unrecognized element found.
             * Will skip over remainder of element later in code.
             */
            annotationFile->addFileReadWarning("Unrecognized element \""
                                               + elementName
                                               + "\" and content ignored.  Updating Workbench "
                                               "may correct this problem.");
            m_stream->skipCurrentElement();
        }
    }
    
    if ( ! annotations.empty()) {
        annotationFile->addAnnotationGroupDuringFileReading(groupType,
                                                            coordSpace,
                                                            tabOrWindowIndex,
                                                            spacerTabIndex,
                                                            mediaFileName,
                                                            histologySpaceKey,
                                                            uniqueKey,
                                                            annotations);
    }
}

/**
 * Read a two dimensional annotation.
 *
 * @param annotationFile
 *     The annotation file being read
 * @param annotationElementName
 *     Name of two-dimensional attribute.
 * @param annotation
 *     Two-dimensional annotation that has its data read.
 */
void
AnnotationFileXmlReader::readOneCoordinateAnnotation(AnnotationFile* annotationFile,
                                                     const QString& annotationElementName,
                                                     AnnotationOneCoordinateShape* annotation)
{
    CaretAssert(annotation);
    
    const QXmlStreamAttributes attributes = m_stream->attributes();
    
    readAnnotationAttributes(annotation,
                             annotationElementName,
                             attributes);
    
    /*
     * Shape width
     */
    annotation->setWidth(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                        annotationElementName,
                                                        ATTRIBUTE_WIDTH));
    /*
     * Shape height
     */
    annotation->setHeight(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                        annotationElementName,
                                                        ATTRIBUTE_HEIGHT));
    /*
     * Shape rotation angle
     */
    annotation->setRotationAngle(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                        annotationElementName,
                                                        ATTRIBUTE_ROTATION_ANGLE));
    
    bool done(false);
    while ( ! done) {
        const QXmlStreamReader::TokenType tokenType(m_stream->readNext());
        if (m_stream->atEnd()) {
            done = true;
        }
        else {
            const QString elementName(m_stream->name().toString());
            switch (tokenType) {
                case QXmlStreamReader::StartElement:
                    if (elementName == ELEMENT_COORDINATE_ONE) {
                        /*
                         * Read the coordinate
                         */
                        const bool readStartElementFlag(false);
                        readCoordinate(ELEMENT_COORDINATE_ONE,
                                       annotation->getCoordinate(),
                                       annotation->getCoordinateSpace(),
                                       readStartElementFlag);
                        m_stream->skipCurrentElement();
                    }
                    else if (elementName == ELEMENT_IMAGE_RGBA_BYTES_IN_BASE64) {
                        AnnotationImage* imageAnn = dynamic_cast<AnnotationImage*>(annotation);
                        if (imageAnn != NULL) {
                            readImageDataElement(imageAnn);
                        }
                        else {
                            /*
                             * Issue warning
                             */
                            annotationFile->addFileReadWarning("Annotation Type \""
                                                               + AnnotationTypeEnum::toGuiName(annotation->getType())
                                                               + " should not contain child element "
                                                               + ELEMENT_IMAGE_RGBA_BYTES_IN_BASE64);
                        }
                        m_stream->skipCurrentElement();
                    }
                    else if (elementName == ELEMENT_TEXT_DATA) {
                        AnnotationText* textAnn = dynamic_cast<AnnotationText*>(annotation);
                        if (textAnn != NULL) {
                            readTextDataElement(textAnn,
                                                annotationElementName);
                        }
                        else {
                            /*
                             * Issue warning
                             */
                            annotationFile->addFileReadWarning("Annotation Type \""
                                                               + AnnotationTypeEnum::toGuiName(annotation->getType())
                                                               + " should not contain child element "
                                                               + ELEMENT_TEXT_DATA);
                        }
                    }
                    else if (elementName == GiftiXmlElements::TAG_METADATA) {
                        /*
                         * Reads through the metadata end element so
                         * DO NOT skip to next element
                         */
                        m_streamHelper->readMetaData(annotation->getMetaData());
                    }
                    else {
                        /*
                         * Issue warning (instead of fatal error) if unrecognized element found.
                         * Will skip over remainder of element later in code.
                         */
                        annotationFile->addFileReadWarning("Unrecognized element \""
                                                           + elementName
                                                           + "\" and content ignored.  Updating Workbench "
                                                           "may correct this problem.");
                        m_stream->skipCurrentElement();
                    }
                    break;
                case QXmlStreamReader::EndElement:
                    if (elementName == annotationElementName) {
                        done = true;
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

/**
 * Read the image annotation element.
 *
 * @param imageAnnotation
 *     Image annotation that has its element data read.
 * @throw
 *     DataFileException
 */
void
AnnotationFileXmlReader::readImageDataElement(AnnotationImage* imageAnnotation)
{
    CaretAssert(imageAnnotation);
    
    const QXmlStreamAttributes attributes = m_stream->attributes();
    
    const int32_t imageWidth = m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                            ELEMENT_IMAGE_RGBA_BYTES_IN_BASE64,
                                                                            ATTRIBUTE_IMAGE_WIDTH);
    const int32_t imageHeight = m_streamHelper->getRequiredAttributeIntValue(attributes,
                                                                            ELEMENT_IMAGE_RGBA_BYTES_IN_BASE64,
                                                                            ATTRIBUTE_IMAGE_HEIGHT);
    const int32_t numberOfBytes = imageWidth * imageHeight * 4;
    
    
    /*
     * Read the image bytes in base64 encoding which will also finish reading through
     * the closing element.
     */
    const QString imageChars = m_stream->readElementText(QXmlStreamReader::ErrorOnUnexpectedElement);
    if (m_stream->hasError()) {
        m_streamHelper->throwDataFileException("There was an error reading the image annotation's image bytes: "
                                               + m_stream->errorString());
    }
    
    QByteArray imageBytes = QByteArray::fromBase64(imageChars.toLatin1());
    if (imageBytes.size() == numberOfBytes) {
        const uint8_t* imageBytesPointer = (const uint8_t*)(imageBytes.data());
        imageAnnotation->setImageBytesRGBA(imageBytesPointer,
                                           imageWidth,
                                           imageHeight);
    }
    else {
        m_streamHelper->throwDataFileException("There was an error reading the image annotations image bytes.  "
                                               "The number of bytes read was " + AString::number(imageBytes.size())
                                               + " but the number of bytes expected was " + AString::number(numberOfBytes));
    }
}

/**
 * Read the text annotation element.
 *
 * @param textAnnotation
 *     Text annotation that has its element data read.
 * @param annotationTextElementName
 *     Name of the annotation element.
 * @throw
 *     DataFileException
 */
void
AnnotationFileXmlReader::readTextDataElement(AnnotationText *textAnnotation,
                                             const QString& annotationTextElementName)
{
    CaretAssert(textAnnotation);
    
    const QXmlStreamAttributes attributes = m_stream->attributes();
    
    bool haveTextColorFlag = false;
    
    {
        /*
         * Background color
         */
        const QString valueString = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                    annotationTextElementName,
                                                                                    ATTRIBUTE_TEXT_CARET_COLOR,
                                                                                    "");
        if ( ! valueString.isEmpty()) {
            bool valid = false;
            CaretColorEnum::Enum value = CaretColorEnum::fromName(valueString,
                                                                  &valid);
            if (valid) {
                textAnnotation->setTextColor(value);
                haveTextColorFlag = true;
            }
            else {
                m_streamHelper->throwDataFileException("Invalid value "
                                                       + valueString
                                                       + " for attribute "
                                                       + ATTRIBUTE_TEXT_CARET_COLOR);
            }
        }
    }
    
    bool haveCustomTextColorFlag = false;
    {
        /*
         * Background custom color
         */
        const QString valueString = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                    annotationTextElementName,
                                                                                    ATTRIBUTE_TEXT_CUSTOM_RGBA,
                                                                                    "");
        if ( ! valueString.isEmpty()) {
            std::vector<float> rgba;
            AString::toNumbers(valueString, rgba);
            if (rgba.size() == 4) {
                textAnnotation->setCustomTextColor(&rgba[0]);
                haveCustomTextColorFlag = true;
            }
            else {
                m_streamHelper->throwDataFileException(ATTRIBUTE_TEXT_CUSTOM_RGBA
                                                       + " must contain 4 elements but "
                                                       + valueString
                                                       + " contains "
                                                       + QString::number(rgba.size())
                                                       + " elements");
            }
        }
    }
    
    if (haveTextColorFlag
        && haveCustomTextColorFlag) {
        /* nothing */
    }
    else {
        /*
         * Older (pre Workbench 1.2) annotations did not have a text color
         * and the text was drawn using the foreground color.
         * So, copy the foreground color to the text color and set
         * the foreground color to none.
         */
        textAnnotation->setTextColor(textAnnotation->getLineColor());
        float rgba[4];
        textAnnotation->getCustomLineColor(rgba);
        textAnnotation->setCustomTextColor(rgba);
        textAnnotation->setLineColor(CaretColorEnum::NONE);
    }
    
    textAnnotation->setBoldStyleEnabled(m_streamHelper->getRequiredAttributeBoolValue(attributes,
                                                                 ELEMENT_TEXT_DATA,
                                                                 ATTRIBUTE_TEXT_FONT_BOLD));
    textAnnotation->setItalicStyleEnabled(m_streamHelper->getRequiredAttributeBoolValue(attributes,
                                                                   ELEMENT_TEXT_DATA,
                                                                   ATTRIBUTE_TEXT_FONT_ITALIC));
    textAnnotation->setUnderlineStyleEnabled(m_streamHelper->getRequiredAttributeBoolValue(attributes,
                                                                   ELEMENT_TEXT_DATA,
                                                                   ATTRIBUTE_TEXT_FONT_UNDERLINE));
    
    {
        const QString defaultValue = AnnotationTextConnectTypeEnum::toName(AnnotationTextConnectTypeEnum::ANNOTATION_TEXT_CONNECT_NONE);
        const QString valueString = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                    ELEMENT_TEXT_DATA,
                                                                                    ATTRIBUTE_TEXT_CONNECT_BRAINORDINATE,
                                                                                    defaultValue);
        bool valid = false;
        AnnotationTextConnectTypeEnum::Enum connectValue = AnnotationTextConnectTypeEnum::fromName(valueString,
                                                                                               &valid);
        if (valid) {
            textAnnotation->setConnectToBrainordinate(connectValue);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_TEXT_CONNECT_BRAINORDINATE);
        }
    }
    
    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                                    ELEMENT_TEXT_DATA,
                                                                                    ATTRIBUTE_TEXT_FONT_NAME);
        bool valid = false;
        AnnotationTextFontNameEnum::Enum fontName = AnnotationTextFontNameEnum::fromName(valueString,
                                                                                 &valid);
        if (valid) {
            textAnnotation->setFont(fontName);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                                   + valueString
                                                   + " for attribute "
                                                   + ATTRIBUTE_TEXT_FONT_NAME);
        }
    }
    
    if ((annotationTextElementName != ELEMENT_PERCENT_SIZE_TEXT)
        && (annotationTextElementName != ELEMENT_PERCENT_WIDTH_SIZE_TEXT)) {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                                    ELEMENT_TEXT_DATA,
                                                                                    ATTRIBUTE_TEXT_FONT_POINT_SIZE,
                                                                                    "fontSize");
        bool valid = false;
        AnnotationTextFontPointSizeEnum::Enum fontPointSize = AnnotationTextFontPointSizeEnum::fromName(valueString,
                                                                                 &valid);
        if (valid) {
            textAnnotation->setFontPointSizeProtected(fontPointSize);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_TEXT_FONT_POINT_SIZE);
        }
    }

    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    ELEMENT_TEXT_DATA,
                                                                    ATTRIBUTE_TEXT_HORIZONTAL_ALIGNMENT);
        bool valid = false;
        AnnotationTextAlignHorizontalEnum::Enum alignment = AnnotationTextAlignHorizontalEnum::fromName(valueString,
                                                                                 &valid);
        if (valid) {
            textAnnotation->setHorizontalAlignment(alignment);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_TEXT_HORIZONTAL_ALIGNMENT);
        }
    }
    
    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    ELEMENT_TEXT_DATA,
                                                                    ATTRIBUTE_TEXT_VERTICAL_ALIGNMENT);
        bool valid = false;
        AnnotationTextAlignVerticalEnum::Enum alignment = AnnotationTextAlignVerticalEnum::fromName(valueString,
                                                                                                        &valid);
        if (valid) {
            textAnnotation->setVerticalAlignment(alignment);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_TEXT_VERTICAL_ALIGNMENT);
        }
    }
    
    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                    ELEMENT_TEXT_DATA,
                                                                    ATTRIBUTE_TEXT_ORIENTATION);
        bool valid = false;
        AnnotationTextOrientationEnum::Enum orientation = AnnotationTextOrientationEnum::fromName(valueString,
                                                                                                    &valid);
        if (valid) {
            textAnnotation->setOrientation(orientation);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                   + valueString
                                   + " for attribute "
                                   + ATTRIBUTE_TEXT_ORIENTATION);
        }
    }
    
    if (annotationTextElementName == ELEMENT_PERCENT_SIZE_TEXT) {
        textAnnotation->setFontPercentViewportSizeProtected(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                                                                  ELEMENT_TEXT_DATA,
                                                                                                  ATTRIBUTE_TEXT_FONT_PERCENT_VIEWPORT_SIZE));
    }
    else if (annotationTextElementName == ELEMENT_PERCENT_WIDTH_SIZE_TEXT) {
        textAnnotation->setFontPercentViewportSizeProtected(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                                                                           ELEMENT_TEXT_DATA,
                                                                                                           ATTRIBUTE_TEXT_FONT_PERCENT_VIEWPORT_SIZE));
    }
    else {
        textAnnotation->setFontPercentViewportSizeProtected(m_streamHelper->getOptionalAttributeFloatValue(attributes,
                                                                                                  ELEMENT_TEXT_DATA,
                                                                                                  ATTRIBUTE_TEXT_FONT_PERCENT_VIEWPORT_SIZE,
                                                                                                  5.0));
    }
    
    
    if (annotationTextElementName == ELEMENT_TEXT_OBSOLETE) {
    }
    else if (annotationTextElementName == ELEMENT_PERCENT_SIZE_TEXT) {
        /*
         * Will cause warning or assertion failure if invalid value
         */
        textAnnotation->setFontPercentViewportSizeProtected(textAnnotation->getFontPercentViewportSizeProtected());
    }
    else if (annotationTextElementName == ELEMENT_PERCENT_WIDTH_SIZE_TEXT) {
        /*
         * Will cause warning or assertion failure if invalid value
         */
        textAnnotation->setFontPercentViewportSizeProtected(textAnnotation->getFontPercentViewportSizeProtected());
    }
    else if (annotationTextElementName == ELEMENT_POINT_SIZE_TEXT) {
        
    }
    else {
        m_streamHelper->throwDataFileException("Unrecognized text element name \""
                                               + annotationTextElementName
                                               + "\" expected "
                                               + ELEMENT_PERCENT_SIZE_TEXT
                                               + " or "
                                               + ELEMENT_PERCENT_WIDTH_SIZE_TEXT
                                               + " or "
                                               + ELEMENT_POINT_SIZE_TEXT);
    }
    
    /*
     * Read the annotation's text which will also finish reading through
     * the closing element.
     */
    const QString textChars = m_stream->readElementText(QXmlStreamReader::ErrorOnUnexpectedElement);
    if (m_stream->hasError()) {
        m_streamHelper->throwDataFileException("There was an error reading the text annotation's characters: "
                               + m_stream->errorString());
    }
    textAnnotation->setText(textChars);
    
    if (m_stream->hasError()) {
        m_streamHelper->throwDataFileException("There was an error reading the annotation file in XML format (reported by QXmlStreamReader): "
                               + m_stream->errorString());
    }
}

/**
 * Read the font attributes from the given XML stream attributes
 * @param fontAttributes
 *    The font attributes
 * @param attributes
 *    The XML stream attributes
 */
void
AnnotationFileXmlReader::readFontAttibutes(AnnotationFontAttributesInterface* fontAttributes,
                                           const AString& elementName,
                                           const QXmlStreamAttributes& attributes)
{
    CaretAssert(fontAttributes);
    
    {
        const QString valueString = m_streamHelper->getRequiredAttributeStringValue(attributes,
                                                                                    elementName,
                                                                                    ATTRIBUTE_TEXT_FONT_NAME);
        bool valid = false;
        AnnotationTextFontNameEnum::Enum fontName = AnnotationTextFontNameEnum::fromName(valueString,
                                                                                         &valid);
        if (valid) {
            fontAttributes->setFont(fontName);
        }
        else {
            m_streamHelper->throwDataFileException("Invalid value "
                                                   + valueString
                                                   + " for attribute "
                                                   + ATTRIBUTE_TEXT_FONT_NAME);
        }
    }

    fontAttributes->setFontPercentViewportSize(m_streamHelper->getRequiredAttributeFloatValue(attributes,
                                                                                              elementName,
                                                                                              ATTRIBUTE_TEXT_FONT_PERCENT_VIEWPORT_SIZE));

    {
        /*
         * Text color
         */
        const QString valueString = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                    elementName,
                                                                                    ATTRIBUTE_TEXT_CARET_COLOR,
                                                                                    "");
        if ( ! valueString.isEmpty()) {
            bool valid = false;
            CaretColorEnum::Enum value = CaretColorEnum::fromName(valueString,
                                                                  &valid);
            if (valid) {
                fontAttributes->setTextColor(value);
            }
            else {
                m_streamHelper->throwDataFileException("Invalid value "
                                                       + valueString
                                                       + " for attribute "
                                                       + ATTRIBUTE_TEXT_CARET_COLOR);
            }
        }
    }

    {
        /*
         * Custom color
         */
        const QString valueString = m_streamHelper->getOptionalAttributeStringValue(attributes,
                                                                                    elementName,
                                                                                    ATTRIBUTE_TEXT_CUSTOM_RGBA,
                                                                                    "");
        if ( ! valueString.isEmpty()) {
            std::vector<float> rgba;
            AString::toNumbers(valueString, rgba);
            if (rgba.size() == 4) {
                fontAttributes->setCustomTextColor(&rgba[0]);
            }
            else {
                m_streamHelper->throwDataFileException(ATTRIBUTE_TEXT_CUSTOM_RGBA
                                                       + " must contain 4 elements but "
                                                       + valueString
                                                       + " contains "
                                                       + QString::number(rgba.size())
                                                       + " elements");
            }
        }
    }

    fontAttributes->setBoldStyleEnabled(m_streamHelper->getRequiredAttributeBoolValue(attributes,
                                                                                      elementName,
                                                                                      ATTRIBUTE_TEXT_FONT_BOLD));
    fontAttributes->setItalicStyleEnabled(m_streamHelper->getRequiredAttributeBoolValue(attributes,
                                                                                        elementName,
                                                                                        ATTRIBUTE_TEXT_FONT_ITALIC));
    fontAttributes->setUnderlineStyleEnabled(m_streamHelper->getRequiredAttributeBoolValue(attributes,
                                                                                           elementName,
                                                                                           ATTRIBUTE_TEXT_FONT_UNDERLINE));
}

