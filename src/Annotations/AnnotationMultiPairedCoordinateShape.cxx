
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

#define __ANNOTATION_MULTI_PAIRED_COORDINATE_SHAPE_DECLARE__
#include "AnnotationMultiPairedCoordinateShape.h"
#undef __ANNOTATION_MULTI_PAIRED_COORDINATE_SHAPE_DECLARE__

#include <algorithm>
#include <cmath>

#include "AnnotationCoordinate.h"
#include "AnnotationPolyhedron.h"
#include "AnnotationSpatialModification.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "MathFunctions.h"
#include "SceneClass.h"
#include "SceneClassArray.h"
#include "SceneClassAssistant.h"

using namespace caret;


    
/**
 * \class caret::AnnotationMultiPairedCoordinateShape
 * \brief Class for annotations that pairs of coordinates
 * \ingroup Annotations
 */

/**
 * Constructor.
 * @param type
 *    Type of annotation.
 * @param attributeDefaultType
 *    Type for attribute defaults
 */
AnnotationMultiPairedCoordinateShape::AnnotationMultiPairedCoordinateShape(const AnnotationTypeEnum::Enum type,
                                                             const AnnotationAttributesDefaultTypeEnum::Enum attributeDefaultType)
: Annotation(type,
             attributeDefaultType)
{
    initializeMembersAnnotationMultiPairedCoordinateShape();
}

/**
 * Destructor.
 */
AnnotationMultiPairedCoordinateShape::~AnnotationMultiPairedCoordinateShape()
{
}

/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
AnnotationMultiPairedCoordinateShape::AnnotationMultiPairedCoordinateShape(const AnnotationMultiPairedCoordinateShape& obj)
: Annotation(obj)
{
    initializeMembersAnnotationMultiPairedCoordinateShape();
    
    this->copyHelperAnnotationMultiPairedCoordinateShape(obj);
}

/**
 * Assignment operator.
 * @param obj
 *    Data copied from obj to this.
 * @return 
 *    Reference to this object.
 */
AnnotationMultiPairedCoordinateShape&
AnnotationMultiPairedCoordinateShape::operator=(const AnnotationMultiPairedCoordinateShape& obj)
{
    if (this != &obj) {
        Annotation::operator=(obj);
        this->copyHelperAnnotationMultiPairedCoordinateShape(obj);
    }
    return *this;    
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void 
AnnotationMultiPairedCoordinateShape::copyHelperAnnotationMultiPairedCoordinateShape(const AnnotationMultiPairedCoordinateShape& obj)
{
    m_coordinates.clear();
    
    for (const auto& ptr : obj.m_coordinates) {
        std::unique_ptr<AnnotationCoordinate> ac(new AnnotationCoordinate(*ptr));
        m_coordinates.push_back(std::move(ac));
    }
    
    setModified();
}

/**
 * Initialize members of this class.
 */
void
AnnotationMultiPairedCoordinateShape::initializeMembersAnnotationMultiPairedCoordinateShape()
{
    m_sceneAssistant.reset(new SceneClassAssistant());
    if (testProperty(Property::SCENE_CONTAINS_ATTRIBUTES)) {
    }
}

/**
 * @return 'this' as a one-dimensional shape. NULL if this is not a one-dimensional shape.
 */
AnnotationMultiPairedCoordinateShape*
AnnotationMultiPairedCoordinateShape::castToMultiPairedCoordinateShape()
{
    return this;
}

/**
 * @return 'this' as a one-dimensional shape. NULL if this is not a one-dimensional shape.
 */
const AnnotationMultiPairedCoordinateShape*
AnnotationMultiPairedCoordinateShape::castToMultiPairedCoordinateShape() const
{
    return this;
}

/**
 * Add a coordinate to this multi coordinate shape
 * @param coord
 *    Coordinate that is added.
 */
void
AnnotationMultiPairedCoordinateShape::addCoordinate(AnnotationCoordinate* coord)
{
    CaretAssert(coord);
    std::unique_ptr<AnnotationCoordinate> ptr(coord);
    m_coordinates.push_back(std::move(ptr));
    setModified();
}

/**
 * Add a coordinate pair to this multi-coordinate shape
 * @param acOne
 *     FIrst coordinate added at the 'middle'
 * @param acTwo
 *     Second coordinate added at the end
 */
void
AnnotationMultiPairedCoordinateShape::addCoordinatePair(AnnotationCoordinate* acOne,
                                                        AnnotationCoordinate* acTwo)
{
    CaretAssert(acOne);
    CaretAssert(acTwo);
    std::unique_ptr<AnnotationCoordinate> ptrOne(acOne);
    std::unique_ptr<AnnotationCoordinate> ptrTwo(acTwo);
    const int32_t numCoords(getNumberOfCoordinates());
    if (numCoords > 0) {
        CaretAssert(MathFunctions::isEvenNumber(numCoords));
        const int32_t halfOffset(numCoords / 2);
        m_coordinates.insert(m_coordinates.begin() + halfOffset,
                             std::move(ptrOne));
        m_coordinates.push_back(std::move(ptrTwo));
    }
    else {
        m_coordinates.push_back(std::move(ptrOne));
        m_coordinates.push_back(std::move(ptrTwo));
    }
    
    CaretAssert(MathFunctions::isEvenNumber(getNumberOfCoordinates()));
}

/**
 * @return Number of coordinates in this annotation
 */
int32_t
AnnotationMultiPairedCoordinateShape::getNumberOfCoordinates() const
{
    return m_coordinates.size();
}

/**
 * @return Coordinate at the given index
 * @param index
 *    Inde of the coordinate
 */
AnnotationCoordinate*
AnnotationMultiPairedCoordinateShape::getCoordinate(const int32_t index)
{
    CaretAssertVectorIndex(m_coordinates, index);
    return m_coordinates[index].get();
}

/**
 * @return Coordinate at the given index
 * @param index
 *    Inde of the coordinate
 */
const AnnotationCoordinate*
AnnotationMultiPairedCoordinateShape::getCoordinate(const int32_t index) const
{
    CaretAssertVectorIndex(m_coordinates, index);
    return m_coordinates[index].get();
}

/**
 * Replace all coordinates in this annotation with copies of the given coordinates
 * @param coordinates
 *    Coordinates that are copied into this annotation
 */
void
AnnotationMultiPairedCoordinateShape::replaceAllCoordinates(const std::vector<std::unique_ptr<const AnnotationCoordinate>>& coordinates)
{
    m_coordinates.clear();
    
    for (const auto& coord : coordinates) {
        AnnotationCoordinate* ac = new AnnotationCoordinate(*coord);
        addCoordinate(ac);
    }
}

/**
 * Insert one coordinate after the first index that is between the two coordinates
 * by the normalized distance
 * @param indexOne
 *    Index of first coordinate
 * @param indexTwo
 *    Index of second coordinate
 * @param normalizedDistanceToNextCoordinate
 *    Normalized distance from first to second coordinate
 */
void
AnnotationMultiPairedCoordinateShape::insertSingleCoordinate(const int32_t indexOne,
                                                             const int32_t indexTwo,
                                                             const float normalizedDistanceToNextCoordinate)
{
    CaretAssertVectorIndex(m_coordinates, indexOne);
    CaretAssertVectorIndex(m_coordinates, indexTwo);
    const Vector3D xyzOne(m_coordinates[indexOne]->getXYZ());
    const Vector3D xyzTwo(m_coordinates[indexTwo]->getXYZ());
    const Vector3D vector((xyzTwo - xyzOne).normal());
    const float distance((xyzTwo - xyzOne).length() * normalizedDistanceToNextCoordinate);
    const Vector3D newXYZ(xyzOne + (vector * distance));
    
    std::unique_ptr<AnnotationCoordinate> newCoord(new AnnotationCoordinate(m_attributeDefaultType));
    newCoord->setXYZ(newXYZ);
    
    const int32_t numCoords(static_cast<int32_t>(m_coordinates.size()));
    const int32_t insertIndex(indexOne + 1);
    if (insertIndex > numCoords) {
        CaretAssertMessage(0, ("Invalid insertion index="
                               + AString::number(insertIndex)
                               + ", number of coords="
                               + AString::number(numCoords)));
        return;
    }
    else if (insertIndex == numCoords) {
        m_coordinates.push_back(std::move(newCoord));
    }
    else {
        m_coordinates.insert(m_coordinates.begin() + insertIndex,
                             std::move(newCoord));
    }
}

/**
 * Insert a new coordinate after the given index. New coordinate is at midpoint to next coordinate.
 * @param insertAfterCoordinateIndex
 *    Insert a coordinate after this coordinate index.
 * @param normalizedDistanceToNextCoordinate
 *    Normalized distance to next coordinate for insertion
 * @param surfaceSpaceVertexIndex
 *    Vertex index for inserting a vertex in a surface space annotation
 */
void
AnnotationMultiPairedCoordinateShape::insertCoordinate(const int32_t insertAfterCoordinateIndex,
                                                       const int32_t /*surfaceSpaceVertexIndex*/,
                                                       const float normalizedDistanceToNextCoordinate)
{
    bool validFlag(false);
    switch (getCoordinateSpace()) {
        case AnnotationCoordinateSpaceEnum::CHART:
            break;
        case AnnotationCoordinateSpaceEnum::HISTOLOGY:
            break;
        case AnnotationCoordinateSpaceEnum::MEDIA_FILE_NAME_AND_PIXEL:
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            validFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            break;
    }
    
    if ( ! validFlag) {
        CaretLogSevere("Multi-paired coord space invalid: "
                       + AnnotationCoordinateSpaceEnum::toName(getCoordinateSpace()));
        return;
    }

    AnnotationPolyhedron* polyhedron = castToPolyhedron();
    if (polyhedron == NULL) {
        AString msg("Shape is not polyhedron.  Has new multi-paired coordinate shape been added?");
        CaretAssertMessage(0, msg);
        CaretLogSevere(msg);
        return;
    }
    
    validFlag = false;
    
    /*
     * Coordinates are in pairs (first set followed by second set)
     */
    const int32_t fullNumCoords = static_cast<int32_t>(m_coordinates.size());
    const int32_t halfNumCoords(fullNumCoords / 2);
    if (halfNumCoords < 2) {
        CaretLogSevere("Multi paired coordinate Shape has invalid number of coordinates="
                       + AString::number(fullNumCoords)
                       + "cannot insert new coordinates.");
        return;
    }
    
    int32_t polyOneIndexOne(-1);
    int32_t polyOneIndexTwo(-1);
    int32_t polyTwoIndexOne(-1);
    int32_t polyTwoIndexTwo(-1);
    if (polyhedron != NULL) {
        /*
         * The polyhedron consists of a sequence of two polygons
         * so two coordinates must be inserted
         */
        if (insertAfterCoordinateIndex < halfNumCoords) {
            /*
             * FIRST POLYGON
             */
            if ((insertAfterCoordinateIndex >= 0)
                && (insertAfterCoordinateIndex < (halfNumCoords - 1))) {
                polyOneIndexOne = insertAfterCoordinateIndex;
                polyOneIndexTwo = insertAfterCoordinateIndex + 1;
                validFlag = true;
            }
            else if (insertAfterCoordinateIndex == -1) {
                polyOneIndexOne = halfNumCoords - 1;
                polyOneIndexTwo = 0;
                validFlag = true;
            }
            else if (insertAfterCoordinateIndex == (halfNumCoords - 1)) {
                polyOneIndexOne = halfNumCoords - 1;
                polyOneIndexTwo = 0;
                validFlag = true;
            }
            polyTwoIndexOne = polyOneIndexOne + halfNumCoords;
            polyTwoIndexTwo = polyOneIndexTwo + halfNumCoords;
        }
        else {
            /*
             * SECOND POLYGON
             */
            if ((insertAfterCoordinateIndex >= halfNumCoords)
                && (insertAfterCoordinateIndex < (fullNumCoords - 1))) {
                polyTwoIndexOne = insertAfterCoordinateIndex;
                polyTwoIndexTwo = insertAfterCoordinateIndex + 1;
                validFlag = true;
            }
            else if (insertAfterCoordinateIndex == -1) {
                polyTwoIndexOne = fullNumCoords - 1;
                polyTwoIndexTwo = halfNumCoords;
                validFlag = true;
            }
            else if (insertAfterCoordinateIndex == (fullNumCoords - 1)) {
                polyTwoIndexOne = fullNumCoords - 1;
                polyTwoIndexTwo = halfNumCoords;
                validFlag = true;
            }
            polyOneIndexOne = polyTwoIndexOne - halfNumCoords;
            polyOneIndexTwo = polyTwoIndexTwo - halfNumCoords;
        }
    }
    else {
        CaretAssert(0);
    }

    if ( ! validFlag) {
        CaretLogSevere("Attempting to insert coordinate after invalid index="
                       + AString::number(insertAfterCoordinateIndex)
                       + " into annotation with coordinate count="
                       + AString::number(fullNumCoords));
        return;
    }
    
    CaretAssertVectorIndex(m_coordinates, polyOneIndexOne);
    CaretAssertVectorIndex(m_coordinates, polyOneIndexTwo);
    CaretAssertVectorIndex(m_coordinates, polyTwoIndexOne);
    CaretAssertVectorIndex(m_coordinates, polyTwoIndexTwo);

    /*
     * Insert coordinate into second polygon
     */
    insertSingleCoordinate(polyTwoIndexOne,
                           polyTwoIndexTwo,
                           normalizedDistanceToNextCoordinate);
    
    /*
     * Insert into first polygon AFTER insertion into first polygon
     * since this insertion shifts the second polygon's coordinates
     */
    insertSingleCoordinate(polyOneIndexOne,
                           polyOneIndexTwo,
                           normalizedDistanceToNextCoordinate);
}

/**
 * Get the coordinates indices for coordinates clockwise and counter-clockwise to the given coordinates
 * @param coordinateIndex
 *    The coordinate index
 * @param clockwiseIndexOut
 *    Output with clockwise index (negative if not valid)
 * @param counterClockwiseIndexOut
 *    Output with counter-clockwise index (negative if not valid)
 * @return True if successful, else false
 */
bool
AnnotationMultiPairedCoordinateShape::getClockwiseAndCounterClockwiseCoordinates(const int32_t coordinateIndex,
                                                                           int32_t& clockwiseIndexOut,
                                                                           int32_t& counterClockwiseIndexOut) const
{
    clockwiseIndexOut        = -1;
    counterClockwiseIndexOut = -1;
    
    bool validSpaceFlag(true);
    switch (getCoordinateSpace()) {
        case AnnotationCoordinateSpaceEnum::CHART:
            break;
        case AnnotationCoordinateSpaceEnum::HISTOLOGY:
            break;
        case AnnotationCoordinateSpaceEnum::MEDIA_FILE_NAME_AND_PIXEL:
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            validSpaceFlag = false;
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            break;
    }
    if ( ! validSpaceFlag) {
        return false;
    }
    const int32_t numCoords(m_coordinates.size());
    if (numCoords < 3) {
        return false;
    }
    
    int32_t i1(0);
    int32_t i2(1);
    int32_t i3(2);
    
    if (numCoords > 3) {
        int32_t index(coordinateIndex);
        if (index == 0) {
            index++;
        }
        const int32_t lastIndex(numCoords - 1);
        if (index >= lastIndex) {
            index = lastIndex - 1;
        }
        i1 = index - 1;
        i2 = index;
        i3 = index + 1;
    }
    
    CaretAssertVectorIndex(m_coordinates, i1);
    CaretAssertVectorIndex(m_coordinates, i2);
    CaretAssertVectorIndex(m_coordinates, i3);

    float normalVector[3];
    if (MathFunctions::normalVector(m_coordinates[i1]->getXYZ(),
                                    m_coordinates[i2]->getXYZ(),
                                    m_coordinates[i3]->getXYZ(),
                                    normalVector)) {
        if (normalVector[2] < 0.0f) {
            std::swap(i1, i3);
        }
        
        clockwiseIndexOut        = i1;
        counterClockwiseIndexOut = i3;
        
        std::cout << "index=" << coordinateIndex
        << " clockwise=" << clockwiseIndexOut
        << " counter-clockwise=" << counterClockwiseIndexOut
        << std::endl << std::flush;
    }
    
    return false;
}

/**
 * Remove the coordinate at the given index
 * @param index
 *    Index of coordinate for removal
 */
void
AnnotationMultiPairedCoordinateShape::removeCoordinateAtIndex(const int32_t index)
{
    const int32_t totalNumCoords(getNumberOfCoordinates());
    const int32_t halfNumCoords(totalNumCoords / 2);
    
    if (index < halfNumCoords) {
        const int32_t indexTwo = index + halfNumCoords;
        CaretAssertVectorIndex(m_coordinates, indexTwo);
        m_coordinates.erase(m_coordinates.begin() + indexTwo);
        CaretAssertVectorIndex(m_coordinates, index);
        m_coordinates.erase(m_coordinates.begin() + index);
    }
    else {
        const int32_t indexTwo = index - halfNumCoords;
        CaretAssertVectorIndex(m_coordinates, index);
        m_coordinates.erase(m_coordinates.begin() + index);
        CaretAssertVectorIndex(m_coordinates, indexTwo);
        m_coordinates.erase(m_coordinates.begin() + indexTwo);
    }
    setModified();
}

/**
 * Remove the coordinate at the given index.  This method should ONLY BE CALLED
 * by UserInputModeAnnotations to erase a coordinate from an annotation that is being drawn
 * BEFORE the second face of the polyhedron is added.
 * @param index
 *    Index of coordinate for removal
 * @param removePairFlag
 *    If true, remove a pair of coordinates.  This should only be true if the annotation is drawn by
 *    adding pairs of coordinates.
 */
void
AnnotationMultiPairedCoordinateShape::removeCoordinateAtIndexByUserInputModeAnnotations(const int32_t index,
                                                                                        const bool removePairFlag)
{
    if (removePairFlag) {
        removeCoordinateAtIndex(index);
    }
    else {
        CaretAssertVectorIndex(m_coordinates, index);
        m_coordinates.erase(m_coordinates.begin() + index);
        setModified();
    }
}


/**
 * @return The surface offset vector type for this annotation.
 */
AnnotationSurfaceOffsetVectorTypeEnum::Enum
AnnotationMultiPairedCoordinateShape::getSurfaceOffsetVectorType() const
{
    AnnotationSurfaceOffsetVectorTypeEnum::Enum offsetType(AnnotationSurfaceOffsetVectorTypeEnum::SURFACE_NORMAL);
    if ( ! m_coordinates.empty()) {
        CaretAssertVectorIndex(m_coordinates, 0);
        offsetType = m_coordinates[0]->getSurfaceOffsetVectorType();
    }
    
    return offsetType;
}

/**
 * Is the object modified?
 * @return true if modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::isModified() const
{
    if (Annotation::isModified()) {
        return true;
    }
    
    for (const auto& ptr : m_coordinates) {
        if (ptr->isModified()) {
            return true;
        }
    }
    
    return false;
}

/**
 * Set the status to unmodified.
 */
void
AnnotationMultiPairedCoordinateShape::clearModified()
{
    Annotation::clearModified();
    
    for (auto& ptr : m_coordinates) {
        ptr->clearModified();
    }
}

/**
 * Apply the coordinates, size, and rotation from the given annotation
 * to this annotation.
 * 
 * @param otherAnnotationIn
 *     The other annotation from which attributes are obtained.
 */
void
AnnotationMultiPairedCoordinateShape::applyCoordinatesSizeAndRotationFromOther(const Annotation* otherAnnotationIn)
{
    CaretAssert(otherAnnotationIn);
    const AnnotationMultiPairedCoordinateShape* otherAnn = otherAnnotationIn->castToMultiPairedCoordinateShape();
    CaretAssert(otherAnn);

    const int32_t otherNumCoords = otherAnn->getNumberOfCoordinates();
    const int32_t myNumCoords    = getNumberOfCoordinates();
    if (myNumCoords == otherNumCoords) {
        for (int32_t i = 0; i < myNumCoords; i++) {
            *getCoordinate(i) = *otherAnn->getCoordinate(i);
        }
    }
    else {
        m_coordinates.clear();
        for (int32_t i = 0; i < otherNumCoords; i++) {
            const AnnotationCoordinate* otherCoord = otherAnn->getCoordinate(i);
            
            std::unique_ptr<AnnotationCoordinate> ac(new AnnotationCoordinate(*otherCoord));
            m_coordinates.push_back(std::move(ac));
        }
    }

    setCoordinateSpace(otherAnn->getCoordinateSpace());
    setTabIndex(otherAnn->getTabIndex());
    setWindowIndex(otherAnn->getWindowIndex());
    
    AnnotationPolyhedron* thisPolyhedron(castToPolyhedron());
    const AnnotationPolyhedron* otherPolyhedron(otherAnn->castToPolyhedron());
    if ((thisPolyhedron != NULL)
        && (otherPolyhedron != NULL)) {
        thisPolyhedron->setPlaneOneNameStereotaxicXYZ(otherPolyhedron->getPlaneOneNameStereotaxicXYZ());
        thisPolyhedron->setPlaneTwoNameStereotaxicXYZ(otherPolyhedron->getPlaneTwoNameStereotaxicXYZ());
    }
    
    setModified();
}

/**
 * Is the given sizing handle valid for this annotation?
 *
 * @sizingHandle
 *    The sizing handle.
 * @return
 *    True if sizing handle valid, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::isSizeHandleValid(const AnnotationSizingHandleTypeEnum::Enum sizingHandle) const
{
    bool xyPlaneFlag = false;
    
    switch (getCoordinateSpace()) {
        case AnnotationCoordinateSpaceEnum::CHART:
            xyPlaneFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::HISTOLOGY:
            xyPlaneFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::MEDIA_FILE_NAME_AND_PIXEL:
            xyPlaneFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            xyPlaneFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            xyPlaneFlag = true;
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            xyPlaneFlag = true;
            break;
    }
    
    bool validFlag = false;
    
    switch (sizingHandle) {
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
            if (xyPlaneFlag) {
                validFlag = true;
            }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
            validFlag = true;
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
            validFlag = true;
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
            validFlag = true;
            break;
    }
    
    return validFlag;
}

/**
 * Apply a spatial modification to an annotation in surface space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationSurfaceSpace(const AnnotationSpatialModification& spatialModification)
{
    bool validFlag = false;
    
    switch (spatialModification.m_sizingHandleType) { 
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
        {
            StructureEnum::Enum structure = StructureEnum::INVALID;
            int32_t surfaceNumberOfNodes  = -1;
            int32_t surfaceNodeIndex      = -1;
            const int32_t coordIndex(spatialModification.m_polyLineCoordinateIndex);
            if ((coordIndex >= 0)
                && (coordIndex < getNumberOfCoordinates())) {
                
                AnnotationCoordinate* coord = getCoordinate(coordIndex);
                CaretAssert(coord);
                coord->getSurfaceSpace(structure,
                                       surfaceNumberOfNodes,
                                       surfaceNodeIndex);
                if (spatialModification.m_surfaceCoordinateAtMouseXY.m_surfaceNodeValid) {
                    if ((spatialModification.m_surfaceCoordinateAtMouseXY.m_surfaceStructure == structure)
                        && (spatialModification.m_surfaceCoordinateAtMouseXY.m_surfaceNumberOfNodes == surfaceNumberOfNodes)) {
                        coord->setSurfaceSpace(spatialModification.m_surfaceCoordinateAtMouseXY.m_surfaceStructure,
                                               spatialModification.m_surfaceCoordinateAtMouseXY.m_surfaceNumberOfNodes,
                                               spatialModification.m_surfaceCoordinateAtMouseXY.m_surfaceNodeIndex);
                        validFlag = true;
                    }
                }
            }
        }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
            break;
    }
    
    
    if (validFlag) {
        setModified();
    }
    
    return validFlag;
}

/**
 * Apply a spatial modification to an annotation in spacer tab space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationSpacerTabSpace(const AnnotationSpatialModification& spatialModification)
{
    return applySpatialModificationTabOrWindowSpace(spatialModification);
}


/**
 * Apply a spatial modification to an annotation in tab or window space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationTabOrWindowSpace(const AnnotationSpatialModification& spatialModification)
{
    bool validFlag = false;
    const int32_t numCoords(getNumberOfCoordinates());
        /*
         * Get all XYZ coordinates
         */
        std::vector<float> allXYZ;
        allXYZ.reserve(numCoords * 3);
        for (int32_t i = 0; i < numCoords; i++) {
            const AnnotationCoordinate* ac(getCoordinate(i));
            float xyz[3];
            ac->getXYZ(xyz);
            float viewportXYZ[3] = { 0.0, 0.0, 0.0 };
            relativeXYZToViewportXYZ(xyz, spatialModification.m_viewportWidth, spatialModification.m_viewportHeight, viewportXYZ);
            allXYZ.push_back(viewportXYZ[0]);
            allXYZ.push_back(viewportXYZ[1]);
            allXYZ.push_back(viewportXYZ[2]);
        }
        
        const float spaceDX = spatialModification.m_mouseDX;
        const float spaceDY = spatialModification.m_mouseDY;
        
        int32_t startIndex(-1);
        int32_t endIndex(-1);
        switch (spatialModification.m_sizingHandleType) {
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
                /*
                 * Moving entire shape (all coordinates change)
                 */
                startIndex = 0;
                endIndex   = numCoords - 1;
                validFlag = true;
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
                /*
                 * Moving one coordinate in the shape
                 */
                if ((spatialModification.m_polyLineCoordinateIndex >= 0)
                    && (spatialModification.m_polyLineCoordinateIndex < numCoords)) {
                    startIndex = spatialModification.m_polyLineCoordinateIndex;
                    endIndex   = spatialModification.m_polyLineCoordinateIndex;
                    validFlag = true;
                }
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
                break;
            case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
                break;
        }
        
        if ((validFlag)
            && (startIndex >= 0)
            && (endIndex >= 0)) {
            
            /*
             * Need to convert back to percentage coords
             */
            for (int32_t i = 0; i < numCoords; i++) {
                const int32_t i3(i * 3);
                float x(allXYZ[i3]);
                float y(allXYZ[i3+1]);
                
                /*
                 * Adjust selected coord(s)
                 */
                if ((i >= startIndex)
                    && (i <= endIndex)) {
                    x += spaceDX;
                    y += spaceDY;
                }
                
                const float percentagNewX = 100.0 * (x / spatialModification.m_viewportWidth);
                const float percentagNewY = 100.0 * (y / spatialModification.m_viewportHeight);
                if ((percentagNewX >= 0.0)
                    && (percentagNewX <= 100.0)
                    && (percentagNewY >= 0.0)
                    && (percentagNewY <= 100.0)) {
                    allXYZ[i3]   = percentagNewX;
                    allXYZ[i3+1] = percentagNewY;
                }
                else {
                    validFlag = false;
                    break;
                }
            }
            
            if (validFlag) {
                /*
                 * Update the coordinates
                 */
                for (int32_t i = 0; i < numCoords; i++) {
                    AnnotationCoordinate* ac(getCoordinate(i));
                    const int32_t i3(i * 3);
                    ac->setXYZ(&allXYZ[i3]);
                }
                if (startIndex == endIndex) {
                    if (spatialModification.isMultiPairedMove()) {
                        movePairedCoordinate(startIndex);
                    }
                }
            }
        }

    if (validFlag) {
        setModified();
    }
    return validFlag;
}

/**
 * Apply a spatial modification to an annotation in chart space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationChartSpace(const AnnotationSpatialModification& spatialModification)
{
    bool validFlag = false;
    const int32_t numCoords(getNumberOfCoordinates());
    int32_t startIndex(-1);
    int32_t endIndex(-1);
    switch (spatialModification.m_sizingHandleType) {
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
            /*
             * Moving entire shape (all coordinates change)
             */
            startIndex = 0;
            endIndex   = numCoords - 1;
            validFlag = true;
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
            /*
             * Moving one coordinate in the shape
             */
            if ((spatialModification.m_polyLineCoordinateIndex >= 0)
                && (spatialModification.m_polyLineCoordinateIndex < numCoords)) {
                startIndex = spatialModification.m_polyLineCoordinateIndex;
                endIndex   = spatialModification.m_polyLineCoordinateIndex;
                validFlag = true;
            }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
            break;
    }
    if ((validFlag)
        && (startIndex >= 0)
        && (endIndex >= 0)) {
        validFlag = false;
        
        if (spatialModification.m_chartCoordAtMouseXY.m_chartXYZValid
            && spatialModification.m_chartCoordAtPreviousMouseXY.m_chartXYZValid) {
            const float dx = spatialModification.m_chartCoordAtMouseXY.m_chartXYZ[0] - spatialModification.m_chartCoordAtPreviousMouseXY.m_chartXYZ[0];
            const float dy = spatialModification.m_chartCoordAtMouseXY.m_chartXYZ[1] - spatialModification.m_chartCoordAtPreviousMouseXY.m_chartXYZ[1];
            const float dz = spatialModification.m_chartCoordAtMouseXY.m_chartXYZ[2] - spatialModification.m_chartCoordAtPreviousMouseXY.m_chartXYZ[2];
            
            for (int32_t i = startIndex; i <= endIndex; i++) {
                AnnotationCoordinate* ac = getCoordinate(i);
                ac->addToXYZ(dx, dy, dz);
            }
            if (startIndex == endIndex) {
                if (spatialModification.isMultiPairedMove()) {
                    movePairedCoordinate(startIndex);
                }
            }
            validFlag = true;
        }
    }
    
    if (validFlag) {
        setModified();
    }
    
    return validFlag;
}

/**
 * Apply a spatial modification to an annotation in histology space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationHistologySpace(const AnnotationSpatialModification& spatialModification)
{
    bool validFlag = false;
    const int32_t numCoords(getNumberOfCoordinates());
    int32_t startIndex(-1);
    int32_t endIndex(-1);
    switch (spatialModification.m_sizingHandleType) {
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
            /*
             * Moving entire shape (all coordinates change)
             */
            startIndex = 0;
            endIndex   = numCoords - 1;
            validFlag = true;
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
            /*
             * Moving one coordinate in the shape
             */
            if ((spatialModification.m_polyLineCoordinateIndex >= 0)
                && (spatialModification.m_polyLineCoordinateIndex < numCoords)) {
                startIndex = spatialModification.m_polyLineCoordinateIndex;
                endIndex   = spatialModification.m_polyLineCoordinateIndex;
                validFlag = true;
            }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
            break;
    }
    if ((validFlag)
        && (startIndex >= 0)
        && (endIndex >= 0)) {
        validFlag = false;
        
        if (spatialModification.m_histologyCoordAtMouseXY.m_histologyXYZValid &&
            spatialModification.m_histologyCoordAtPreviousMouseXY.m_histologyXYZValid) {
            const float dx = spatialModification.m_histologyCoordAtMouseXY.m_histologyXYZ[0] - spatialModification.m_histologyCoordAtPreviousMouseXY.m_histologyXYZ[0];
            const float dy = spatialModification.m_histologyCoordAtMouseXY.m_histologyXYZ[1] - spatialModification.m_histologyCoordAtPreviousMouseXY.m_histologyXYZ[1];
            const float dz = spatialModification.m_histologyCoordAtMouseXY.m_histologyXYZ[2] - spatialModification.m_histologyCoordAtPreviousMouseXY.m_histologyXYZ[2];
            
            for (int32_t i = startIndex; i <= endIndex; i++) {
                AnnotationCoordinate* ac = getCoordinate(i);
                ac->addToXYZ(dx, dy, dz);
            }
            if (startIndex == endIndex) {
                if (spatialModification.isMultiPairedMove()) {
                    movePairedCoordinate(startIndex);
                }
            }
            validFlag = true;
        }
    }
    
    if (validFlag) {
        setModified();
    }
    
    return validFlag;
}

/**
 * Apply a spatial modification to an annotation in media space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationMediaSpace(const AnnotationSpatialModification& spatialModification)
{
    bool validFlag = false;
    const int32_t numCoords(getNumberOfCoordinates());
    int32_t startIndex(-1);
    int32_t endIndex(-1);
    switch (spatialModification.m_sizingHandleType) {
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
            /*
             * Moving entire shape (all coordinates change)
             */
            startIndex = 0;
            endIndex   = numCoords - 1;
            validFlag = true;
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
            /*
             * Moving one coordinate in the shape
             */
            if ((spatialModification.m_polyLineCoordinateIndex >= 0)
                && (spatialModification.m_polyLineCoordinateIndex < numCoords)) {
                startIndex = spatialModification.m_polyLineCoordinateIndex;
                endIndex   = spatialModification.m_polyLineCoordinateIndex;
                validFlag = true;
            }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
            break;
    }
    if ((validFlag)
        && (startIndex >= 0)
        && (endIndex >= 0)) {
        validFlag = false;
        
        if (spatialModification.m_mediaCoordAtMouseXY.m_mediaXYZValid &&
            spatialModification.m_mediaCoordAtPreviousMouseXY.m_mediaXYZValid) {
            const float dx = spatialModification.m_mediaCoordAtMouseXY.m_mediaXYZ[0] - spatialModification.m_mediaCoordAtPreviousMouseXY.m_mediaXYZ[0];
            const float dy = spatialModification.m_mediaCoordAtMouseXY.m_mediaXYZ[1] - spatialModification.m_mediaCoordAtPreviousMouseXY.m_mediaXYZ[1];
            const float dz = spatialModification.m_mediaCoordAtMouseXY.m_mediaXYZ[2] - spatialModification.m_mediaCoordAtPreviousMouseXY.m_mediaXYZ[2];
            
            for (int32_t i = startIndex; i <= endIndex; i++) {
                AnnotationCoordinate* ac = getCoordinate(i);
                ac->addToXYZ(dx, dy, dz);
            }
            if (startIndex == endIndex) {
                if (spatialModification.isMultiPairedMove()) {
                    movePairedCoordinate(startIndex);
                }
            }
            validFlag = true;
        }
    }
    
    if (validFlag) {
        setModified();
    }
    
    return validFlag;
}

/**
 * Apply a spatial modification to an annotation in stereotaxic space.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModificationStereotaxicSpace(const AnnotationSpatialModification& spatialModification)
{
    bool modifiedFlag(false);
    switch (spatialModification.m_sizingHandleType) {
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
            /* No dragging entire annotation in stereotaxic space */
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_EDITABLE_POLY_LINE_COORDINATE:
        {
            const int32_t coordIndex(spatialModification.m_polyLineCoordinateIndex);
            const int32_t numCoords(getNumberOfCoordinates());
            if ((coordIndex >= 0)
                && (coordIndex < numCoords)) {
                AnnotationCoordinate* ac = getCoordinate(coordIndex);
                ac->setXYZ(spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ);
                if (spatialModification.isMultiPairedMove()) {
                    movePairedCoordinate(coordIndex);
                }
                setModified();
                modifiedFlag = true;
            }
        }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NOT_EDITABLE_POLY_LINE_COORDINATE:
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_ONE:
        {
            AnnotationPolyhedron* polyhedron(castToPolyhedron());
            if (polyhedron != NULL) {
                const Vector3D xyz(spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ[0],
                                   spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ[1],
                                   spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ[2]);
                polyhedron->setPlaneOneNameStereotaxicXYZ(xyz);
                modifiedFlag = true;
            }
        }
            break;
        case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLYHEDRON_TEXT_COORDINATE_TWO:
        {
            AnnotationPolyhedron* polyhedron(castToPolyhedron());
            if (polyhedron != NULL) {
                const Vector3D xyz(spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ[0],
                                   spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ[1],
                                   spatialModification.m_stereotaxicCoordinateAtMouseXY.m_stereotaxicXYZ[2]);
                polyhedron->setPlaneTwoNameStereotaxicXYZ(xyz);
                modifiedFlag = true;
            }
        }
            break;
    }

    return modifiedFlag;
}

/**
 * Move the paired coordinate that is paired to the given coordinate
 * @param coordIndex
 *    Index of the coordinate
 */
void
AnnotationMultiPairedCoordinateShape::movePairedCoordinate(const int32_t coordIndex)
{
    AnnotationPolyhedron* polyhedron = castToPolyhedron();
    if (polyhedron == NULL) {
        AString msg("Shape is not polyhedron.  Has new multi-paired coordinate shape been added?");
        CaretAssertMessage(0, msg);
        CaretLogSevere(msg);
        return;
    }

    /*
     * Coords are in PAIRS
     */
    Plane plane;
    const int32_t numCoords(getNumberOfCoordinates());
    const int32_t numCoordPairs(numCoords / 2);
    int32_t pairedCoordIndex(-1);
    if ((coordIndex >= 0)
        && (coordIndex < numCoordPairs)) {
        pairedCoordIndex = coordIndex + numCoordPairs;
        plane = polyhedron->getPlaneTwo();
    }
    else if ((coordIndex >= numCoordPairs)
             && (coordIndex < numCoords)) {
        pairedCoordIndex = coordIndex - numCoordPairs;
        plane = polyhedron->getPlaneOne();
    }
    else {
        CaretLogSevere("Invalid coordinate index="
                       + AString::number(coordIndex));
    }
    
    if ((pairedCoordIndex >= 0)
        && (pairedCoordIndex < numCoords)) {
        if (plane.isValidPlane()) {
            const AnnotationCoordinate* ac(getCoordinate(coordIndex));
            Vector3D xyz(ac->getXYZ());
            const Vector3D xyzTwo(plane.projectPointToPlane(xyz));
            AnnotationCoordinate* acPair(getCoordinate(pairedCoordIndex));
            acPair->setXYZ(xyzTwo);
        }
        else {
            CaretLogSevere("Plane is invalid for moving paired coordinate");
        }
    }
}

/**
 * Apply a spatial modification to an annotation.
 *
 * @param spatialModification
 *     Contains information about the spatial modification.
 * @return
 *     True if the annotation was modified, else false.
 */
bool
AnnotationMultiPairedCoordinateShape::applySpatialModification(const AnnotationSpatialModification& spatialModification)
{
    if ( ! isSizeHandleValid(spatialModification.m_sizingHandleType)) {
        return false;
    }
    
    switch (getCoordinateSpace()) {
        case AnnotationCoordinateSpaceEnum::CHART:
            return applySpatialModificationChartSpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::HISTOLOGY:
            return applySpatialModificationHistologySpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::MEDIA_FILE_NAME_AND_PIXEL:
            return applySpatialModificationMediaSpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            return applySpatialModificationSpacerTabSpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            return applySpatialModificationStereotaxicSpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            return applySpatialModificationSurfaceSpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            return applySpatialModificationTabOrWindowSpace(spatialModification);
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            return applySpatialModificationTabOrWindowSpace(spatialModification);
            break;
    }
    
    return false;
}

/**
 * Save subclass data to the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass to which data members should be added.  Will always
 *     be valid (non-NULL).
 */
void
AnnotationMultiPairedCoordinateShape::saveSubClassDataToScene(const SceneAttributes* sceneAttributes,
                                            SceneClass* sceneClass)
{
    m_sceneAssistant->saveMembers(sceneAttributes,
                                  sceneClass);
}

/**
 * Restore file data from the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass for the instance of a class that implements
 *     this interface.  Will NEVER be NULL.
 */
void
AnnotationMultiPairedCoordinateShape::restoreSubClassDataFromScene(const SceneAttributes* sceneAttributes,
                                                 const SceneClass* sceneClass)
{
    m_sceneAssistant->restoreMembers(sceneAttributes,
                                     sceneClass);
}

