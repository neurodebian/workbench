
/*LICENSE_START*/
/*
 *  Copyright (C) 2016 Washington University School of Medicine
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

#define __CHART_TWO_COMPOUND_DATA_TYPE_DECLARE__
#include "ChartTwoCompoundDataType.h"
#undef __CHART_TWO_COMPOUND_DATA_TYPE_DECLARE__

#include "CaretAssert.h"
#include "SceneClass.h"
#include "SceneClassAssistant.h"

using namespace caret;


    
/**
 * \class caret::ChartTwoCompoundDataType 
 * \brief Compound type for type of chart
 * \ingroup Charting
 *
 * A compound type for matching chart types.
 * Equality is determined using several pieces of information.
 */

/**
 * Constructor.  It is best to use the "new instance" methods.
 *
 * @param chartDataType
 *    Chart data type.
 * @param histogramNumberOfBuckets
 *     Number of bins in the histogram.
 * @param lineChartUnitsAxisX
 *    Line chart x-axis units.
 * @param lineChartUnitsAxisY
 *    Line chart y-axis units.
 * @param lineChartNumberOfElementsAxisX
 *    Line chart x-axis number of elements.
 * @param matrixNumberOfRows
 *    Matrix number of rows.
 * @param matrixNumberOfColumns 
 *    Matrix number of columns.
 */
ChartTwoCompoundDataType::ChartTwoCompoundDataType(const ChartTwoDataTypeEnum::Enum chartDataType,
                                                   const int32_t histogramNumberOfBuckets,
                                                   const CaretUnitsTypeEnum::Enum lineChartUnitsAxisX,
                                                   const CaretUnitsTypeEnum::Enum lineChartUnitsAxisY,
                                                   const int32_t lineChartNumberOfElementsAxisX,
                                                   const int32_t matrixNumberOfRows,
                                                   const int32_t matrixNumberOfColumns)
: CaretObject()
{
    initializeChartTwoCompoundDataType();
    
    m_chartDataType                  = chartDataType;
    m_histogramNumberOfBuckets       = histogramNumberOfBuckets;
    m_lineChartUnitsAxisX            = lineChartUnitsAxisX;
    m_lineChartUnitsAxisY            = lineChartUnitsAxisY;
    m_lineChartNumberOfElementsAxisX = lineChartNumberOfElementsAxisX;
    m_matrixNumberOfRows             = matrixNumberOfRows;
    m_matrixNumberOfColumns          = matrixNumberOfColumns;
}

/**
 * Default constructor of invalid instance.
 */
ChartTwoCompoundDataType::ChartTwoCompoundDataType()
{
    initializeChartTwoCompoundDataType();
}

/**
 * Destructor.
 */
ChartTwoCompoundDataType::~ChartTwoCompoundDataType()
{
    delete m_sceneAssistant;
}

/**
 * @return A new instance for a histogram chart.
 *
 * @param histogramNumberOfBins
 *     Number of bins in the histogram.
 */
ChartTwoCompoundDataType
ChartTwoCompoundDataType::newInstanceForHistogram(const int32_t histogramNumberOfBins)
{
    return ChartTwoCompoundDataType(ChartTwoDataTypeEnum::CHART_DATA_TYPE_HISTOGRAM,
                                    histogramNumberOfBins,
                                    CaretUnitsTypeEnum::NONE,
                                    CaretUnitsTypeEnum::NONE,
                                    0,
                                    0,
                                    0);
}

/**
 * @return A new instance for a line-layer chart.
 *
 * @param lineChartUnitsAxisX
 *    Line chart x-axis units.
 * @param lineChartUnitsAxisY
 *    Units for the Y-axis
 * @param lineChartNumberOfElementsAxisX
 *    Line chart x-axis number of elements.
 */
ChartTwoCompoundDataType
ChartTwoCompoundDataType::newInstanceForLineLayer(const CaretUnitsTypeEnum::Enum lineChartUnitsAxisX,
                                                  const CaretUnitsTypeEnum::Enum lineChartUnitsAxisY,
                                                  const int32_t lineChartNumberOfElementsAxisX)
{
    return ChartTwoCompoundDataType(ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_LAYER,
                                    0,
                                    lineChartUnitsAxisX,
                                    lineChartUnitsAxisY,
                                    lineChartNumberOfElementsAxisX,
                                    0,
                                    0);
}

/**
 * @return A new instance for a line-series chart.
 *
 * @param lineChartUnitsAxisX
 *    Line chart x-axis units.
 * @param lineChartUnitsAxisY
 *    Units for the Y-axis
 * @param lineChartNumberOfElementsAxisX
 *    Line chart x-axis number of elements.
 */
ChartTwoCompoundDataType
ChartTwoCompoundDataType::newInstanceForLineSeries(const CaretUnitsTypeEnum::Enum lineChartUnitsAxisX,
                                                   const CaretUnitsTypeEnum::Enum lineChartUnitsAxisY,
                                                   const int32_t lineChartNumberOfElementsAxisX)
{
    return ChartTwoCompoundDataType(ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_SERIES,
                                    0,
                                    lineChartUnitsAxisX,
                                    lineChartUnitsAxisY,
                                    lineChartNumberOfElementsAxisX,
                                    0,
                                    0);
}

/**
 * @return A new instance for a matrix chart.
 *
 * @param lineChartUnitsAxisX
 *    Units for the X-axis
 * @param lineChartUnitsAxisY
 *    Units for the Y-axis
 * @param matrixNumberOfRows
 *    Matrix number of rows.
 * @param matrixNumberOfColumns
 *    Matrix number of columns.
 */
ChartTwoCompoundDataType
ChartTwoCompoundDataType::newInstanceForMatrix(const CaretUnitsTypeEnum::Enum lineChartUnitsAxisX,
                                               const CaretUnitsTypeEnum::Enum lineChartUnitsAxisY,
                                               const int32_t matrixNumberOfRows,
                                               const int32_t matrixNumberOfColumns)
{
    return ChartTwoCompoundDataType(ChartTwoDataTypeEnum::CHART_DATA_TYPE_MATRIX,
                                    0,
                                    lineChartUnitsAxisX,
                                    lineChartUnitsAxisY,
                                    0,
                                    matrixNumberOfRows,
                                    matrixNumberOfColumns);
}


/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
ChartTwoCompoundDataType::ChartTwoCompoundDataType(const ChartTwoCompoundDataType& obj)
: CaretObject(obj),
SceneableInterface(obj)
{
    initializeChartTwoCompoundDataType();
    this->copyHelperChartTwoCompoundDataType(obj);
}

/**
 * Assignment operator.
 * @param obj
 *    Data copied from obj to this.
 * @return 
 *    Reference to this object.
 */
ChartTwoCompoundDataType&
ChartTwoCompoundDataType::operator=(const ChartTwoCompoundDataType& obj)
{
    if (this != &obj) {
        CaretObject::operator=(obj);
        this->copyHelperChartTwoCompoundDataType(obj);
    }
    return *this;    
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void 
ChartTwoCompoundDataType::copyHelperChartTwoCompoundDataType(const ChartTwoCompoundDataType& obj)
{
    m_chartDataType = obj.m_chartDataType;
    m_lineChartUnitsAxisX = obj.m_lineChartUnitsAxisX;
    m_lineChartNumberOfElementsAxisX = obj.m_lineChartNumberOfElementsAxisX;
    m_matrixNumberOfRows = obj.m_matrixNumberOfRows;
    m_matrixNumberOfColumns = obj.m_matrixNumberOfColumns;
}

/**
 * Initialize an invalid instance of this class.
 */
void
ChartTwoCompoundDataType::initializeChartTwoCompoundDataType()
{
    m_chartDataType                  = ChartTwoDataTypeEnum::CHART_DATA_TYPE_INVALID;
    m_lineChartUnitsAxisX            = CaretUnitsTypeEnum::NONE;
    m_lineChartUnitsAxisY            = CaretUnitsTypeEnum::NONE;
    m_lineChartNumberOfElementsAxisX = 0;
    m_matrixNumberOfRows             = 0;
    m_matrixNumberOfColumns          = 0;
    
    m_sceneAssistant = new SceneClassAssistant();
    m_sceneAssistant->add<ChartTwoDataTypeEnum, ChartTwoDataTypeEnum::Enum>("m_chartDataType",
                                                                            &m_chartDataType);
    m_sceneAssistant->add<CaretUnitsTypeEnum, CaretUnitsTypeEnum::Enum>("m_lineChartUnitsAxisX",
                                                                        &m_lineChartUnitsAxisX);
    m_sceneAssistant->add<CaretUnitsTypeEnum, CaretUnitsTypeEnum::Enum>("m_lineChartUnitsAxisY",
                                                                        &m_lineChartUnitsAxisY);
    m_sceneAssistant->add("m_lineChartNumberOfElementsAxisX",
                          &m_lineChartNumberOfElementsAxisX);
    m_sceneAssistant->add("m_matrixNumberOfRows",
                          &m_matrixNumberOfRows);
    m_sceneAssistant->add("m_matrixNumberOfColumns",
                          &m_matrixNumberOfColumns);
}


/**
 * Equality operator.
 *
 * Charts must be same type.
 * Line charts must be same units OR same number of elements in X-axis (Y-axis ignored)
 * Matrix charts must have same number of rows and columns.
 *
 * @param obj
 *    Instance compared to this for equality.
 * @return 
 *    True if this instance and 'obj' instance are considered equal.
 */
bool
ChartTwoCompoundDataType::operator==(const ChartTwoCompoundDataType& obj) const
{
    if (this == &obj) {
        return true;    
    }
    
    if (m_chartDataType == obj.m_chartDataType) {
        switch (m_chartDataType) {
            case ChartTwoDataTypeEnum::CHART_DATA_TYPE_INVALID:
                return true;
                break;
            case ChartTwoDataTypeEnum::CHART_DATA_TYPE_HISTOGRAM:
                return true;
                break;
            case ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_LAYER:
            {
                if (m_lineChartUnitsAxisX == obj.m_lineChartUnitsAxisX) {
                    return true;
                }
                if (m_lineChartNumberOfElementsAxisX == obj.m_lineChartNumberOfElementsAxisX) {
                    return true;
                }
            }
                break;
            case ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_SERIES:
            {
                if (m_lineChartUnitsAxisX == obj.m_lineChartUnitsAxisX) {
                    return true;
                }
                if (m_lineChartNumberOfElementsAxisX == obj.m_lineChartNumberOfElementsAxisX) {
                    return true;
                }
            }
                break;
            case ChartTwoDataTypeEnum::CHART_DATA_TYPE_MATRIX:
                if ((m_matrixNumberOfRows == obj.m_matrixNumberOfRows)
                    && (m_matrixNumberOfColumns == obj.m_matrixNumberOfColumns)) {
                    return true;
                }
                break;
        }
    }

    return false;    
}

/**
 * @return Type of chart data.
 */
ChartTwoDataTypeEnum::Enum
ChartTwoCompoundDataType::getChartTwoDataType() const
{
    return m_chartDataType;
}

/**
 * @return Histogram number of buckets.
 */
int32_t
ChartTwoCompoundDataType::getHistogramNumberOfBuckets() const
{
    return m_histogramNumberOfBuckets;
}

/**
 * @return Line chart X-axis units.
 */
CaretUnitsTypeEnum::Enum
ChartTwoCompoundDataType::getLineChartUnitsAxisX() const
{
    return m_lineChartUnitsAxisX;
}

/**
 * @return Line chart Y-axis units.
 */
CaretUnitsTypeEnum::Enum
ChartTwoCompoundDataType::getLineChartUnitsAxisY() const
{
    return m_lineChartUnitsAxisY;
}

/**
 * @return Line chart number of element in the X-axis.
 */
int32_t
ChartTwoCompoundDataType::getLineChartNumberOfElementsAxisX() const
{
    return m_lineChartNumberOfElementsAxisX;
}

/**
 * @return Matrix number of rows.
 */
int32_t
ChartTwoCompoundDataType::getMatrixNumberOfRows() const
{
    return m_matrixNumberOfRows;
}

/**
 * @return Matrix number of columns.
 */
int32_t
ChartTwoCompoundDataType::getMatrixNumberOfColumns() const
{
    return m_matrixNumberOfColumns;
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
ChartTwoCompoundDataType::toString() const
{
    const AString indent = "    ";
    AString text = ("type="
                    + ChartTwoDataTypeEnum::toGuiName(m_chartDataType));
    
    switch (m_chartDataType) {
        case ChartTwoDataTypeEnum::CHART_DATA_TYPE_INVALID:
            break;
        case ChartTwoDataTypeEnum::CHART_DATA_TYPE_HISTOGRAM:
            break;
        case ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_LAYER:
            text.appendWithNewLine(indent
                                   + "x-units="
                                   + CaretUnitsTypeEnum::toName(m_lineChartUnitsAxisX));
            text.appendWithNewLine(indent
                                   + "y-units="
                                   + CaretUnitsTypeEnum::toName(m_lineChartUnitsAxisY));
            text.appendWithNewLine(indent
                                   + "x-elements="
                                   + AString::number(m_lineChartNumberOfElementsAxisX));
            break;
        case ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_SERIES:
            text.appendWithNewLine(indent
                                   + "x-units="
                                   + CaretUnitsTypeEnum::toName(m_lineChartUnitsAxisX));
            text.appendWithNewLine(indent
                                   + "y-units="
                                   + CaretUnitsTypeEnum::toName(m_lineChartUnitsAxisY));
            text.appendWithNewLine(indent
                                   + "x-elements="
                                   + AString::number(m_lineChartNumberOfElementsAxisX));
            break;
        case ChartTwoDataTypeEnum::CHART_DATA_TYPE_MATRIX:
            text.appendWithNewLine(indent
                                   + "x-units="
                                   + CaretUnitsTypeEnum::toName(m_lineChartUnitsAxisX));
            text.appendWithNewLine(indent
                                   + "y-units="
                                   + CaretUnitsTypeEnum::toName(m_lineChartUnitsAxisY));
            text.appendWithNewLine(indent
                                   + "rows="
                                   + AString::number(m_matrixNumberOfRows));
            text.appendWithNewLine(indent
                                   + "columns="
                                   + AString::number(m_matrixNumberOfColumns));
            break;
    }
    
    return text;
}

/**
 * Save information specific to this type of model to the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    saving the scene.
 *
 * @param instanceName
 *    Name of instance in the scene.
 */
SceneClass*
ChartTwoCompoundDataType::saveToScene(const SceneAttributes* sceneAttributes,
                                 const AString& instanceName)
{
    SceneClass* sceneClass = new SceneClass(instanceName,
                                            "ChartTwoCompoundDataType",
                                            1);
    m_sceneAssistant->saveMembers(sceneAttributes,
                                  sceneClass);
    
    // Uncomment if sub-classes must save to scene
    //saveSubClassDataToScene(sceneAttributes,
    //                        sceneClass);
    
    return sceneClass;
}

/**
 * Restore information specific to the type of model from the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass from which model specific information is obtained.
 */
void
ChartTwoCompoundDataType::restoreFromScene(const SceneAttributes* sceneAttributes,
                                      const SceneClass* sceneClass)
{
    if (sceneClass == NULL) {
        return;
    }
    
    m_sceneAssistant->restoreMembers(sceneAttributes,
                                     sceneClass);    
    
    //Uncomment if sub-classes must restore from scene
    //restoreSubClassDataFromScene(sceneAttributes,
    //                             sceneClass);
    
}

