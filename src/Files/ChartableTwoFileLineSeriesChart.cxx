
/*LICENSE_START*/
/*
 *  Copyright (C) 2017 Washington University School of Medicine
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

#define __CHARTABLE_TWO_FILE_LINE_SERIES_CHART_DECLARE__
#include "ChartableTwoFileLineSeriesChart.h"
#undef __CHARTABLE_TWO_FILE_LINE_SERIES_CHART_DECLARE__

#include "CaretAssert.h"
#include "ChartTwoDataCartesian.h"
#include "ChartTwoLineSeriesHistory.h"
#include "CiftiMappableDataFile.h"
#include "CiftiScalarDataSeriesFile.h"
#include "EventChartTwoLoadLineSeriesData.h"
#include "EventManager.h"
#include "SceneClass.h"
#include "SceneClassAssistant.h"

using namespace caret;
    
/**
 * \class caret::ChartableTwoFileLineSeriesChart 
 * \brief Implementation of base chart delegate for line series charts.
 * \ingroup Files
 */

/**
 * Constructor.
 *
 * @param lineSeriesContentType
 *     Content type of the line series data.
 * @param parentCaretMappableDataFile
 *     Parent caret mappable data file that this delegate supports.
 */
ChartableTwoFileLineSeriesChart::ChartableTwoFileLineSeriesChart(const ChartTwoLineSeriesContentTypeEnum::Enum lineSeriesContentType,
                                                                                 CaretMappableDataFile* parentCaretMappableDataFile)
: ChartableTwoFileBaseChart(ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_SERIES,
                                    parentCaretMappableDataFile),
m_lineSeriesContentType(lineSeriesContentType)
{
    m_sceneAssistant = std::unique_ptr<SceneClassAssistant>(new SceneClassAssistant());
    
    CaretUnitsTypeEnum::Enum xAxisUnits = CaretUnitsTypeEnum::NONE;
    int32_t xAxisNumberOfElements = 0;
    
    CaretMappableDataFile* cmdf = getCaretMappableDataFile();
    CaretAssert(cmdf);
    
    switch (lineSeriesContentType) {
        case ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_UNSUPPORTED:
            break;
        case ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_BRAINORDINATE_DATA:
            if (cmdf->getNumberOfMaps() > 1) {
                const NiftiTimeUnitsEnum::Enum mapUnits = cmdf->getMapIntervalUnits();
                xAxisUnits = CaretUnitsTypeEnum::NONE;
                switch (mapUnits) {
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_HZ:
                        xAxisUnits = CaretUnitsTypeEnum::HERTZ;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_MSEC:
                        xAxisUnits = CaretUnitsTypeEnum::SECONDS;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_PPM:
                        xAxisUnits = CaretUnitsTypeEnum::PARTS_PER_MILLION;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_SEC:
                        xAxisUnits = CaretUnitsTypeEnum::SECONDS;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_USEC:
                        xAxisUnits = CaretUnitsTypeEnum::SECONDS;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_UNKNOWN:
                        break;
                }
                xAxisNumberOfElements = cmdf->getNumberOfMaps();
            }
            break;
        case ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_ROW_SCALAR_DATA:
        {
            const CiftiMappableDataFile* ciftiMapFile = getCiftiMappableDataFile();
            CaretAssert(ciftiMapFile);
            std::vector<int64_t> dims;
            ciftiMapFile->getMapDimensions(dims);
            CaretAssertVectorIndex(dims, 1);
            const int32_t numCols = dims[0];
            const int32_t numRows = dims[1];
            
            if ((numRows > 0)
                && (numCols > 1)) {
                const NiftiTimeUnitsEnum::Enum mapUnits = ciftiMapFile->getMapIntervalUnits();
                xAxisUnits = CaretUnitsTypeEnum::NONE;
                switch (mapUnits) {
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_HZ:
                        xAxisUnits = CaretUnitsTypeEnum::HERTZ;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_MSEC:
                        xAxisUnits = CaretUnitsTypeEnum::SECONDS;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_PPM:
                        xAxisUnits = CaretUnitsTypeEnum::PARTS_PER_MILLION;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_SEC:
                        xAxisUnits = CaretUnitsTypeEnum::SECONDS;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_USEC:
                        xAxisUnits = CaretUnitsTypeEnum::SECONDS;
                        break;
                    case NiftiTimeUnitsEnum::NIFTI_UNITS_UNKNOWN:
                        break;
                }
                xAxisNumberOfElements = numCols;
            }
        }            break;
    }

    m_lineChartHistory = std::unique_ptr<ChartTwoLineSeriesHistory>(new ChartTwoLineSeriesHistory(GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN));
    m_sceneAssistant->add("m_lineChartHistory",
                           "ChartTwoLineSeriesHistory",
                            m_lineChartHistory.get());

    /*
     * Must have two or more elements
     */
    if (xAxisNumberOfElements <= 1) {
        m_lineSeriesContentType = ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_UNSUPPORTED;
    }
    
    /*
     * For Cifti Files, use units from dimensions
     */
    CaretUnitsTypeEnum::Enum yAxisUnits = CaretUnitsTypeEnum::NONE;
    {
        const CiftiMappableDataFile* ciftiMapFile = getCiftiMappableDataFile();
        if (ciftiMapFile != NULL) {
            float start(0.0), step(0.0);
            ciftiMapFile->getDimensionUnits(CiftiXML::ALONG_ROW,
                                            xAxisUnits,
                                            start,
                                            step);
            ciftiMapFile->getDimensionUnits(CiftiXML::ALONG_COLUMN,
                                            yAxisUnits,
                                            start,
                                            step);
        }
    }
    
    updateChartTwoCompoundDataTypeAfterFileChanges(ChartTwoCompoundDataType::newInstanceForLineSeries(xAxisUnits,
                                                                                                      yAxisUnits,
                                                                                                      xAxisNumberOfElements));
    
    if (m_lineSeriesContentType != ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_UNSUPPORTED) {
        EventManager::get()->addEventListener(this,
                                              EventTypeEnum::EVENT_CHART_TWO_LOAD_LINE_SERIES_DATA);
    }
}

/**
 * Destructor.
 */
ChartableTwoFileLineSeriesChart::~ChartableTwoFileLineSeriesChart()
{
    EventManager::get()->removeEventFromListener(this,
                                                 EventTypeEnum::EVENT_CHART_TWO_LOAD_LINE_SERIES_DATA);
}

/**
 * Receive an event.
 *
 * @param event
 *    An event for which this instance is listening.
 */
void
ChartableTwoFileLineSeriesChart::receiveEvent(Event* event)
{
    if (event->getEventType() == EventTypeEnum::EVENT_CHART_TWO_LOAD_LINE_SERIES_DATA) {
        const EventChartTwoLoadLineSeriesData* lineSeriesDataEvent = dynamic_cast<EventChartTwoLoadLineSeriesData*>(event);
        CaretAssert(lineSeriesDataEvent);
        
        if (m_lineSeriesContentType != ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_UNSUPPORTED) {
            const std::vector<int32_t> tabIndicesForLoading = getTabIndicesForLoadingData(lineSeriesDataEvent->getValidTabIndices());
            if ( ! tabIndicesForLoading.empty()) {
                
                loadLineCharts(lineSeriesDataEvent);
            }
        }
        
        event->setEventProcessed();
    }
    else {
        ChartableTwoFileBaseChart::receiveEvent(event);
    }
}

/**
 * Load line-series charts.
 *
 * @param lineSeriesDataEvent
 *     Event indicating data for loading.
 */
void
ChartableTwoFileLineSeriesChart::loadLineCharts(const EventChartTwoLoadLineSeriesData* lineSeriesDataEvent)
{
    const MapFileDataSelector mapFileDataSelector = lineSeriesDataEvent->getMapFileDataSelector();
    
    int32_t scalarRowIndex = -1;
    bool loadDataFlag = false;
    switch (m_lineSeriesContentType) {
        case ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_UNSUPPORTED:
            return;
            break;
        case ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_BRAINORDINATE_DATA:
            switch (mapFileDataSelector.getDataSelectionType()) {
                case MapFileDataSelector::DataSelectionType::INVALID:
                    break;
                case MapFileDataSelector::DataSelectionType::COLUMN_DATA:
                    break;
                case MapFileDataSelector::DataSelectionType::ROW_DATA:
                    break;
                case MapFileDataSelector::DataSelectionType::SURFACE_VERTEX:
                    loadDataFlag = true;
                    break;
                case MapFileDataSelector::DataSelectionType::SURFACE_VERTICES_AVERAGE:
                    loadDataFlag = true;
                    break;
                case MapFileDataSelector::DataSelectionType::VOLUME_XYZ:
                    loadDataFlag = true;
                    break;
            }
            break;
        case ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_ROW_SCALAR_DATA:
            switch (mapFileDataSelector.getDataSelectionType()) {
                case MapFileDataSelector::DataSelectionType::INVALID:
                    break;
                case MapFileDataSelector::DataSelectionType::COLUMN_DATA:
                {
                    CaretMappableDataFile* mapFile = NULL;
                    AString mapFileName;
                    int32_t columnIndex = -1;
                    mapFileDataSelector.getColumnIndex(mapFile, mapFileName, columnIndex);
                    if (mapFile == getCaretMappableDataFile()) {
                        loadDataFlag = true;
                    }
                }
                    break;
                case MapFileDataSelector::DataSelectionType::ROW_DATA:
                {
                    CaretMappableDataFile* mapFile = NULL;
                    AString mapFileName;
                    int32_t rowIndex = -1;
                    mapFileDataSelector.getRowIndex(mapFile, mapFileName, rowIndex);
                    if (mapFile == getCaretMappableDataFile()) {
                        loadDataFlag = true;
                    }
                    scalarRowIndex = rowIndex;
                }
                    break;
                case MapFileDataSelector::DataSelectionType::SURFACE_VERTEX:
                    break;
                case MapFileDataSelector::DataSelectionType::SURFACE_VERTICES_AVERAGE:
                    break;
                case MapFileDataSelector::DataSelectionType::VOLUME_XYZ:
                    break;
            }
            break;
    }
    
    if (loadDataFlag) {
        std::vector<float> data;
        getCaretMappableDataFile()->getDataForSelector(mapFileDataSelector,
                                                       data);
        if ( ! data.empty()) {
            const CaretUnitsTypeEnum::Enum xUnits = getChartTwoCompoundDataType()->getLineChartUnitsAxisX();
            CaretAssert(getChartTwoCompoundDataType()->getLineChartNumberOfElementsAxisX() == static_cast<int32_t>(data.size()));
            ChartTwoDataCartesian* cartesianData = new ChartTwoDataCartesian(ChartTwoDataTypeEnum::CHART_DATA_TYPE_LINE_SERIES,
                                                                             xUnits,
                                                                             CaretUnitsTypeEnum::NONE,
                                                                             m_lineChartHistory->getDefaultGraphicsPrimitiveType());
            cartesianData->setMapFileDataSelector(mapFileDataSelector);
            
            float x = 0.0f;
            float xStep = 0.0f;
            getCaretMappableDataFile()->getMapIntervalStartAndStep(x, xStep);
            
            const int32_t numData = static_cast<int32_t>(data.size());
            for (int32_t i = 0; i < numData; i++) {
                CaretAssertVectorIndex(data, i);
                cartesianData->addPoint(x, data[i]);
                x += xStep;
            }
            
            m_lineChartHistory->addHistoryItem(cartesianData);
            
            if (scalarRowIndex >= 0) {
                CaretMappableDataFile* mapFile = getCaretMappableDataFile();
                CaretAssert(mapFile);
                CiftiScalarDataSeriesFile* scalarFile = dynamic_cast<CiftiScalarDataSeriesFile*>(mapFile);
                if (scalarFile != NULL) {
                    for (auto tabIndex : lineSeriesDataEvent->getValidTabIndices()) {
                        scalarFile->setSelectedMapIndex(tabIndex, scalarRowIndex);
                    }
                }
            }
        }
    }
}

/**
 * Load data for the given row or column.
 *
 * @param tabIndex
 *     Index of tab.
 * @param rowOrColumnIndex
 *     Index of row/column for loading.
 */
void
ChartableTwoFileLineSeriesChart::loadDataForRowOrColumn(const int32_t tabIndex,
                                                        const int32_t rowOrColumnIndex)
{
    std::vector<int32_t> tabIndices;
    tabIndices.push_back(tabIndex);

    MapFileDataSelector mapFileSelector;
    mapFileSelector.setRowIndex(getCaretMappableDataFile(),
                                "",
                                rowOrColumnIndex);

    EventChartTwoLoadLineSeriesData lineSeriesDataEvent(tabIndices,
                                    mapFileSelector);
    loadLineCharts(&lineSeriesDataEvent);
}

/**
 * Find tab indices for which user has loading of data enabled in the given valid tab indices.
 *
 * @param validTabIndices
 *     All tabs that are valid.
 * @return 
 *     Tab indices for which data should be loaded.
 */
std::vector<int32_t>
ChartableTwoFileLineSeriesChart::getTabIndicesForLoadingData(const std::vector<int32_t>& validTabIndices) const
{
    std::vector<int32_t> tabIndicesOut;
    
    if ( ! validTabIndices.empty()) {
        if (m_lineChartHistory->isLoadingEnabled()) {
            tabIndicesOut = validTabIndices;
        }
    }
    
    return tabIndicesOut;
}


/**
 * @return Content type of the line series data.
 */
ChartTwoLineSeriesContentTypeEnum::Enum
ChartableTwoFileLineSeriesChart::getLineSeriesContentType() const
{
    return m_lineSeriesContentType;
}

/**
 * @return History of line charts.
 */
ChartTwoLineSeriesHistory*
ChartableTwoFileLineSeriesChart::getHistory()
{
    return m_lineChartHistory.get();
}

/**
 * @return History of line charts (const method)
 */
const ChartTwoLineSeriesHistory*
ChartableTwoFileLineSeriesChart::getHistory() const
{
    return m_lineChartHistory.get();
}


/**
 * @return Is this charting valid ?
 */
bool
ChartableTwoFileLineSeriesChart::isValid() const
{
    return (m_lineSeriesContentType != ChartTwoLineSeriesContentTypeEnum::LINE_SERIES_CONTENT_UNSUPPORTED);
}

/**
 * @retrurn Is this charting empty (no data at this time)
 */
bool
ChartableTwoFileLineSeriesChart::isEmpty() const
{
    if ( ! isValid()) {
        return true;
    }
    
    return (m_lineChartHistory->getHistoryCount() <= 0);
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
ChartableTwoFileLineSeriesChart::saveSubClassDataToScene(const SceneAttributes* sceneAttributes,
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
ChartableTwoFileLineSeriesChart::restoreSubClassDataFromScene(const SceneAttributes* sceneAttributes,
                                                 const SceneClass* sceneClass)
{
    m_sceneAssistant->restoreMembers(sceneAttributes,
                                     sceneClass);
}

