#ifndef __MODEL_CHART_TWO_H__
#define __MODEL_CHART_TWO_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2014  Washington University School of Medicine
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

#include <memory>

#include <QSharedPointer>
#include <QWeakPointer>

#include <list>
#include <set>
#include <map>

#include "CaretDataFileSelectionModel.h"
#include "ChartTwoDataTypeEnum.h"
#include "EventListenerInterface.h"
#include "MapYokingGroupEnum.h"
#include "Model.h"
#include "StructureEnum.h"

namespace caret {

    class CaretDataFileSelectionModel;
    class ChartTwoOverlaySetArray;
    class CiftiConnectivityMatrixParcelFile;
    class CiftiMappableDataFile;
    class ModelChart;
    class OverlaySetArray;
    class SurfaceFile;
    
    /// Controls the display of a chart.
    class ModelChartTwo : public Model, public EventListenerInterface  {
        
    public:
        ModelChartTwo(Brain* brain);
        
        virtual ~ModelChartTwo();
        
        virtual void initializeOverlays() override;
        
        virtual AString getNameForGUI(const bool includeStructureFlag) const override;
        
        virtual AString getNameForBrowserTab() const override;
        
        void loadChartDataForSurfaceNode(const StructureEnum::Enum structure,
                                         const int32_t surfaceNumberOfNodes,
                                         const int32_t nodeIndex);
        
        void loadAverageChartDataForSurfaceNodes(const StructureEnum::Enum structure,
                                                 const int32_t surfaceNumberOfNodes,
                                                 const std::vector<int32_t>& nodeIndices);
        
        void loadChartDataForVoxelAtCoordinate(const float xyz[3]);
        
        void loadChartDataForCiftiMappableFileRow(CiftiMappableDataFile* ciftiMapFile,
                                                  const int32_t rowIndex);
        
        void loadChartDataForYokedScalarDataSeriesFiles(const MapYokingGroupEnum::Enum mapYokingGroup,
                                                        const int32_t mapIndex);
        
        void selectRowColumnInYokedScalarDataSeriesFileOverlay(const MapYokingGroupEnum::Enum mapYokingGroup,
                                                               const int32_t mapIndex);
        
        virtual OverlaySet* getOverlaySet(const int tabIndex) override;
        
        virtual const OverlaySet* getOverlaySet(const int tabIndex) const override;
        
        virtual ChartTwoOverlaySet* getChartTwoOverlaySet(const int tabIndex) override;
        
        virtual const ChartTwoOverlaySet* getChartTwoOverlaySet(const int tabIndex) const override;
        
        std::vector<ChartTwoOverlaySet*> getAllChartTwoOverlaySets(const int32_t tabIndex) const;
        
        virtual void receiveEvent(Event* event) override;
        
        void getValidChartTwoDataTypes(std::vector<ChartTwoDataTypeEnum::Enum>& validChartDataTypesOut) const;
        
        ChartTwoDataTypeEnum::Enum getSelectedChartTwoDataType(const int32_t tabIndex) const;
        
        void setSelectedChartTwoDataType(const int32_t tabIndex,
                                         const ChartTwoDataTypeEnum::Enum dataType);
        
        virtual void copyTabContent(const int32_t sourceTabIndex,
                                    const int32_t destinationTabIndex) override;
        
        void reset();
        
        void restoreSceneFromChartOneModel(ModelChart* modelChartOne);
        
        void copyChartTwoCartesianAxes(const int32_t sourceTabIndex,
                                            const int32_t destinationTabIndex);
    
    protected:
        virtual void saveModelSpecificInformationToScene(const SceneAttributes* sceneAttributes,
                                                         SceneClass* sceneClass) override;
        
        virtual void restoreModelSpecificInformationFromScene(const SceneAttributes* sceneAttributes,
                                                              const SceneClass* sceneClass) override;
        
    private:
        ModelChartTwo(const ModelChartTwo&);
        
        ModelChartTwo& operator=(const ModelChartTwo&);
        
        void initializeCharts();
        
        void removeAllCharts();
        
        void updateChartOverlaySets(const int32_t tabIndex);
        
        void restoreMatrixChartFromChartOneModel(ModelChart* modelChartOne,
                                                 const int32_t tabIndex);
        
        void restoreLineSeriesChartFromChartOneModel(ModelChart* modelChartOne);
        
        /** Overlays sets for this model and for each tab */
        OverlaySetArray* m_overlaySetArray;
        
        mutable ChartTwoDataTypeEnum::Enum m_selectedChartTwoDataType[BrainConstants::MAXIMUM_NUMBER_OF_BROWSER_TABS];
        
        /** Chart Overlay sets for histogram data type */
        std::unique_ptr<ChartTwoOverlaySetArray> m_histogramChartOverlaySetArray;
        
        /** Chart Overlay sets for line series data type */
        std::unique_ptr<ChartTwoOverlaySetArray> m_lineSeriesChartOverlaySetArray;
        
        /** Chart Overlay sets for layer data type */
        std::unique_ptr<ChartTwoOverlaySetArray> m_lineLayerChartOverlaySetArray;
        
        /** Chart Overlay sets for matrix data type */
        std::unique_ptr<ChartTwoOverlaySetArray> m_matrixChartOverlaySetArray;
        
        std::unique_ptr<SceneClassAssistant> m_sceneAssistant;
    };

} // namespace

#endif // __MODEL_CHART_TWO_H__
