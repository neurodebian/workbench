
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

#define __IDENTIFICATION_FORMATTED_TEXT_GENERATOR_DECLARE__
#include "IdentificationFormattedTextGenerator.h"
#undef __IDENTIFICATION_FORMATTED_TEXT_GENERATOR_DECLARE__

#include "Border.h"
#include "Brain.h"
#include "BrainStructure.h"
#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "CaretMappableDataFile.h"
#include "ChartDataCartesian.h"
#include "ChartDataSource.h"
#include "ChartModelDataSeries.h"
#include "ChartTwoDataCartesian.h"
#include "ChartableMatrixInterface.h"
#include "ChartableTwoFileDelegate.h"
#include "ChartableTwoFileHistogramChart.h"
#include "ChartableTwoFileLineSeriesChart.h"
#include "ChartableTwoFileMatrixChart.h"
#include "CiftiMappableConnectivityMatrixDataFile.h"
#include "CiftiMappableDataFile.h"
#include "CaretVolumeExtension.h"
#include "DataToolTipsManager.h"
#include "EventBrowserTabGetAll.h"
#include "EventBrowserTabIndicesGetAllViewed.h"
#include "EventCaretMappableDataFilesAndMapsInDisplayedOverlays.h"
#include "EventCaretMappableDataFilesGet.h"
#include "EventManager.h"
#include "FileIdentificationAttributes.h"
#include "FileInformation.h"
#include "FociFile.h"
#include "Focus.h"
#include "GiftiLabel.h"
#include "GraphicsPrimitive.h"
#include "GraphicsPrimitiveV3f.h"
#include "Histogram.h"
#include "HtmlTableBuilder.h"
#include "IdentificationFilter.h"
#include "IdentificationManager.h"
#include "ImageFile.h"
#include "MapFileDataSelector.h"
#include "MetricDynamicConnectivityFile.h"
#include "OverlaySet.h"
#include "SelectionItemBorderSurface.h"
#include "SelectionItemChartDataSeries.h"
#include "SelectionItemChartFrequencySeries.h"
#include "SelectionItemChartMatrix.h"
#include "SelectionItemCiftiConnectivityMatrixRowColumn.h"
#include "SelectionItemChartTimeSeries.h"
#include "SelectionItemChartTwoHistogram.h"
#include "SelectionItemChartTwoLineSeries.h"
#include "SelectionItemChartTwoMatrix.h"
#include "SelectionItemFocusSurface.h"
#include "SelectionItemFocusVolume.h"
#include "SelectionItemImage.h"
#include "SelectionItemSurfaceNode.h"
#include "SelectionItemVoxel.h"
#include "SelectionManager.h"
#include "IdentificationStringBuilder.h"
#include "LabelFile.h"
#include "MetricFile.h"
#include "Surface.h"
#include "VolumeDynamicConnectivityFile.h"
#include "SurfaceProjectedItem.h"
#include "SurfaceProjectionBarycentric.h"
#include "SurfaceProjectionVanEssen.h"
#include "VolumeFile.h"

using namespace caret;


    
/**
 * \class IdentificationFormattedTextGenerator
 * \brief Creates text describing selected data.
 *
 * Examine the selected data and generate descriptive text.
 */

/**
 * Constructor.
 */
IdentificationFormattedTextGenerator::IdentificationFormattedTextGenerator()
: CaretObject()
{
    
}

/**
 * Destructor.
 */
IdentificationFormattedTextGenerator::~IdentificationFormattedTextGenerator()
{
    
}

/**
 * Create identification text from selection in the identification manager.
 * @param idManselectionManagerager
 *    Selection manager containing selection.
 * @param brain
 *    The brain.
 * @param tabIndex
 *    Index of tab where identification took place
 */
AString 
IdentificationFormattedTextGenerator::createIdentificationText(const SelectionManager* selectionManager,
                                                               const Brain* brain,
                                                               const int32_t tabIndex) const
{
    CaretAssert(selectionManager);
    CaretAssert(brain);
    
    const IdentificationManager* idManager = brain->getIdentificationManager();
    const IdentificationFilter* filter = idManager->getIdentificationFilter();
    
    IdentificationStringBuilder idText;
    std::unique_ptr<HtmlTableBuilder> chartHtmlTableBuilder = createHtmlTableBuilder(3);
    chartHtmlTableBuilder->setTitleBold("Charts");
    std::unique_ptr<HtmlTableBuilder> geometryHtmlTableBuilder = createHtmlTableBuilder(3);
    geometryHtmlTableBuilder->setTitleBold("Geometry");
    std::unique_ptr<HtmlTableBuilder> imageHtmlTableBuilder = createHtmlTableBuilder(2);
    imageHtmlTableBuilder->setTitlePlain("Image");
    std::unique_ptr<HtmlTableBuilder> labelHtmlTableBuilder = createHtmlTableBuilder(2);
    labelHtmlTableBuilder->setTitlePlain("Labels");
    std::unique_ptr<HtmlTableBuilder> layersHtmlTableBuilder = createHtmlTableBuilder(3);
    layersHtmlTableBuilder->setTitlePlain("Layers");
    std::unique_ptr<HtmlTableBuilder> scalarHtmlTableBuilder = createHtmlTableBuilder(2);
    scalarHtmlTableBuilder->setTitlePlain("Scalars");

    const SelectionItemSurfaceNode* surfaceID = selectionManager->getSurfaceNodeIdentification();
    
    this->generateSurfaceVertexIdentificationText(*geometryHtmlTableBuilder,
                                                  brain,
                                                  surfaceID);
    this->generateVolumeVoxelIdentificationText(*geometryHtmlTableBuilder,
                                                brain,
                                                selectionManager->getVoxelIdentification());

    const std::vector<EventCaretMappableDataFilesAndMapsInDisplayedOverlays::FileInfo> displayedFiles = getFilesForIdentification(filter,
                                                                                                                                  tabIndex);

    for (auto fileInfo : displayedFiles) {
        switch (fileInfo.m_overlayType) {
            case EventCaretMappableDataFilesAndMapsInDisplayedOverlays::OverlayType::BRAINORDINATE:
                if (surfaceID->isValid()) {
                    if (fileInfo.m_mapFile->isSurfaceMappable()) {
                        this->generateSurfaceDataIdentificationText(*labelHtmlTableBuilder,
                                                                    *scalarHtmlTableBuilder,
                                                                    fileInfo.m_mapFile,
                                                                    fileInfo.m_mapIndices,
                                                                    brain,
                                                                    surfaceID);
                    }
                }
                if (selectionManager->getVoxelIdentification()->isValid()) {
                    if (fileInfo.m_mapFile->isVolumeMappable()) {
                        this->generateVolumeDataIdentificationText(*labelHtmlTableBuilder,
                                                                   *scalarHtmlTableBuilder,
                                                                   fileInfo.m_mapFile,
                                                                   fileInfo.m_mapIndices,
                                                                   brain,
                                                                   selectionManager->getVoxelIdentification());
                    }
                }
                break;
            case EventCaretMappableDataFilesAndMapsInDisplayedOverlays::OverlayType::CHART_ONE:
                break;
            case EventCaretMappableDataFilesAndMapsInDisplayedOverlays::OverlayType::CHART_TWO:
                this->generateChartTwoHistogramIdentificationText(*chartHtmlTableBuilder,
                                                                  idText,
                                                                  selectionManager->getChartTwoHistogramIdentification(),
                                                                  fileInfo.m_mapFile,
                                                                  fileInfo.m_mapIndices,
                                                                  false);
                
                this->generateChartTwoLineSeriesIdentificationText(*chartHtmlTableBuilder,
                                                                   idText,
                                                                   selectionManager->getChartTwoLineSeriesIdentification(),
                                                                   fileInfo.m_mapFile,
                                                                   fileInfo.m_mapIndices,
                                                                   false);
                
                this->generateChartTwoMatrixIdentificationText(*chartHtmlTableBuilder,
                                                               idText,
                                                               selectionManager->getChartTwoMatrixIdentification(),
                                                               fileInfo.m_mapFile,
                                                               fileInfo.m_mapIndices,
                                                               false);
                
                break;
        }
    }
    

    if (filter->isShowFociEnabled()) {
        this->generateSurfaceFocusIdentifcationText(*layersHtmlTableBuilder,
                                                    selectionManager->getSurfaceFocusIdentification(),
                                                    false);
        this->generateVolumeFocusIdentifcationText(*layersHtmlTableBuilder,
                                                   selectionManager->getVolumeFocusIdentification());
    }

    if (filter->isShowBorderEnabled()) {
        this->generateSurfaceBorderIdentifcationText(*layersHtmlTableBuilder,
                                                     idText,
                                                     selectionManager->getSurfaceBorderIdentification(),
                                                     false);
    }
    
    
    this->generateChartDataSeriesIdentificationText(*chartHtmlTableBuilder,
                                                    selectionManager->getChartDataSeriesIdentification());
    
    this->generateChartFrequencySeriesIdentificationText(*chartHtmlTableBuilder,
                                                         selectionManager->getChartFrequencySeriesIdentification());
    
    this->generateChartTimeSeriesIdentificationText(*chartHtmlTableBuilder,
                                                    selectionManager->getChartTimeSeriesIdentification());
    
    this->generateChartMatrixIdentificationText(*chartHtmlTableBuilder,
                                                selectionManager->getChartMatrixIdentification());
    
    this->generateCiftiConnectivityMatrixIdentificationText(*chartHtmlTableBuilder,
                                                            selectionManager->getCiftiConnectivityMatrixRowColumnIdentification());
    
    this->generateImageIdentificationText(*imageHtmlTableBuilder,
                                          selectionManager->getImageIdentification());
    
    AString textOut;
    textOut.append(geometryHtmlTableBuilder->getAsHtmlTable());
    textOut.append(labelHtmlTableBuilder->getAsHtmlTable());
    textOut.append(scalarHtmlTableBuilder->getAsHtmlTable());
    textOut.append(layersHtmlTableBuilder->getAsHtmlTable());
    textOut.append(chartHtmlTableBuilder->getAsHtmlTable());
    textOut.append(imageHtmlTableBuilder->getAsHtmlTable());
    return textOut;
    
//    QString htmlString("<html>\n"
//                       "<head>\n"
//                       "<style type=\"text/css\">\n"
//                       "table { border-width: 1px; border-style: ridge; }\n"
//                       "table td { padding: 3px; }"
//                       "table th { padding: 2px; }"
//                       "</style>\n"
//                       "</head>\n"
//                       "</body>\n"
//                       + htmlTableBuilder->getAsHtmlTable()
//                       + "\n</body></html>");
//    return htmlString;
}

/**
 * @return Files for identification
 * @param filter
 * Identification filter
 * @param tabIndex
 * Index of tab where ID took place
 */
std::vector<EventCaretMappableDataFilesAndMapsInDisplayedOverlays::FileInfo>
IdentificationFormattedTextGenerator::getFilesForIdentification(const IdentificationFilter* filter,
                                                                const int32_t tabIndex) const
{
    /**
     * Event gets files from enabled overlays in the viewed tab(s)
     */
    EventCaretMappableDataFilesAndMapsInDisplayedOverlays overlayFilesEvent;
    switch (filter->getTabFiltering()) {
        case IdentificationFilterTabSelectionEnum::ALL_DISPLAYED_TABS:
        {
            EventBrowserTabIndicesGetAllViewed viewedTabsEvent;
            EventManager::get()->sendEvent(viewedTabsEvent.getPointer());
            const std::vector<int32_t> tabs = viewedTabsEvent.getAllBrowserTabIndices();
            std::set<int32_t> tabsSet(tabs.begin(),
                                      tabs.end());
            overlayFilesEvent.setTabIndicesConstraint(tabsSet);
        }
            break;
        case IdentificationFilterTabSelectionEnum::MOUSE_CLICKED_TAB:
            if (tabIndex >= 0) {
                std::set<int32_t> tabIndices { tabIndex };
                overlayFilesEvent.setTabIndicesConstraint(tabIndices);
            }
            break;
    }
    EventManager::get()->sendEvent(overlayFilesEvent.getPointer());

    /*
     * Event gets all files
     */
    EventCaretMappableDataFilesGet mapFilesEvent;
    EventManager::get()->sendEvent(mapFilesEvent.getPointer());
    std::vector<CaretMappableDataFile*> allMapFiles;
    mapFilesEvent.getAllFiles(allMapFiles);
    
    /*
     * Test all files for any that the users has enabled for identification
     * on the Identification Dialog's Filtering tab and add them to the
     * overlay files event
     */
    for (auto mapFile : allMapFiles) {
        CaretAssert(mapFile);
        const FileIdentificationAttributes* fileAtts = mapFile->getFileIdentificationAttributes();
        CaretAssert(fileAtts);
        if (fileAtts->isEnabled()) {
            switch (fileAtts->getMapSelectionMode()) {
                case FileIdentificationMapSelectionEnum::ALL:
                {
                    const int32_t numMaps = mapFile->getNumberOfMaps();
                    for (int32_t iMap = 0; iMap < numMaps; iMap++) {
                        overlayFilesEvent.addBrainordinateFileAndMap(mapFile,
                                                                     iMap,
                                                                     tabIndex);
                    }
                }
                    break;
                case FileIdentificationMapSelectionEnum::SELECTED:
                    overlayFilesEvent.addBrainordinateFileAndMap(mapFile,
                                                                 fileAtts->getMapIndex(),
                                                                 tabIndex);
                    break;
            }
        }
    }
    
    /*
     * Get the displayed and user selected file and return them
     */
    std::vector<EventCaretMappableDataFilesAndMapsInDisplayedOverlays::FileInfo> displayedFiles = overlayFilesEvent.getFilesAndMaps();
    return displayedFiles;
}

/**
 * Get text for the tooltip for a selected node.
 *
 * @param brain
 *     The Brain.
 * @param browserTab
 *     Browser tab in which tooltip is displayed
 * @param selectionManager
 *     The selection manager.
 * @param dataToolTipsManager
 *     The data tooltips manager
 * @param idText
 *     String builder for identification text.
 */
AString
IdentificationFormattedTextGenerator::createToolTipText(const Brain* brain,
                                               const BrowserTabContent* browserTab,
                                               const SelectionManager* selectionManager,
                                               const DataToolTipsManager* dataToolTipsManager) const
{
    CaretAssert(brain);
    CaretAssert(browserTab);
    CaretAssert(selectionManager);
    CaretAssert(dataToolTipsManager);
    
    const SelectionItemSurfaceNode* selectedNode = selectionManager->getSurfaceNodeIdentification();
    const SelectionItemVoxel* selectedVoxel = selectionManager->getVoxelIdentification();
    
    IdentificationStringBuilder idText;
    
    if (selectedNode->isValid()) {
        generateSurfaceToolTip(brain,
                               browserTab,
                               selectionManager,
                               dataToolTipsManager,
                               idText);
    }
    else if (selectedVoxel->isValid()) {
        generateVolumeToolTip(browserTab,
                              selectionManager,
                              dataToolTipsManager,
                              idText);
    }
    else {
        generateChartToolTip(selectionManager,
                             dataToolTipsManager,
                             idText);
    }

    AString text;
    if (idText.length() > 0) {
       text = idText.toStringWithHtmlBodyForToolTip();
    }

    return text;
}


/**
 * Generate identification text for volume voxel identification.
 *
 * @param htmlTableBuilder
 *     Html table builder for identification text.
 * @param brain
 *     The brain.
 * @param idVolumeVoxel
 *     Information for volume voxel ID.
 */
void
IdentificationFormattedTextGenerator::generateVolumeVoxelIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                              const Brain* /*brain*/,
                                                              const SelectionItemVoxel* idVolumeVoxel) const
{
    if (idVolumeVoxel->isValid() == false) {
        return;
    }
    
    int64_t ijk[3];
    const VolumeMappableInterface* idVolumeFile = idVolumeVoxel->getVolumeFile();
    idVolumeVoxel->getVoxelIJK(ijk);
    float x, y, z;
    idVolumeFile->indexToSpace(ijk[0], ijk[1], ijk[2], x, y, z);
    
    const QString xyzText(AString::number(x, 'f', 2)
                          + ", "
                          + AString::number(y, 'f', 2)
                          + ", "
                          + AString::number(z, 'f', 2));
    
    const QString ijkText("Voxel IJK ("
                          + AString::fromNumbers(ijk, 3, ", ")
                          + ")");

    QString filename;
    const CaretDataFile* caretDataFile = dynamic_cast<const CaretDataFile*>(idVolumeFile);
    if (caretDataFile != NULL) {
        filename = caretDataFile->getFileNameNoPath();
    }
    
    htmlTableBuilder.addHeaderRow(ijkText,
                                  xyzText,
                                  filename);
}

/**
 * Generate identification text for volume data identification.
 *
 * @param labelHtmlTableBuilder
 *     HTML table builder for label identification text.
 * @param scalarHtmlTableBuilder
 *     HTML table builder for scalar identification text.
 * @param mapFile
 *     FIle for generating identification
 * @param mapIndices
 *     Indices of map
 * @param brain
 *     The brain.
 * @param idVolumeVoxel
 *     Information for volume voxel ID.
 */
void
IdentificationFormattedTextGenerator::generateVolumeDataIdentificationText(HtmlTableBuilder& labelHtmlTableBuilder,
                                                                           HtmlTableBuilder& scalarHtmlTableBuilder,
                                                                           CaretMappableDataFile* mapFile,
                                                                           const std::set<int32_t>& mapIndicesSet,
                                                                           const Brain* brain,
                                                                           const SelectionItemVoxel* idVolumeVoxel) const
{
    if (idVolumeVoxel->isValid() == false) {
        return;
    }
    
    int64_t ijk[3];
    const VolumeMappableInterface* idVolumeFile = idVolumeVoxel->getVolumeFile();
    idVolumeVoxel->getVoxelIJK(ijk);
    float x, y, z;
    idVolumeFile->indexToSpace(ijk[0], ijk[1], ijk[2], x, y, z);
    const float xyz[3] = { x, y, z };

    /*
     * Get all volume files
     */
    std::vector<const VolumeMappableInterface*> volumeInterfaces;
    const int32_t numVolumeFiles = brain->getNumberOfVolumeFiles();
    for (int32_t i = 0; i < numVolumeFiles; i++) {
        const VolumeFile* vf = brain->getVolumeFile(i);
        volumeInterfaces.push_back(vf);

        const VolumeDynamicConnectivityFile* volDynConnFile = vf->getVolumeDynamicConnectivityFile();
        if (volDynConnFile != NULL) {
            if (volDynConnFile->isDataValid()) {
                volumeInterfaces.push_back(volDynConnFile);
            }
        }
    }

    /*
     * Get the CIFTI files that are volume mappable
     */
    std::vector<CiftiMappableDataFile*> allCiftiMappableDataFiles;
    brain->getAllCiftiMappableDataFiles(allCiftiMappableDataFiles);
    for (std::vector<CiftiMappableDataFile*>::iterator ciftiMapIter = allCiftiMappableDataFiles.begin();
         ciftiMapIter != allCiftiMappableDataFiles.end();
         ciftiMapIter++) {
        const CiftiMappableDataFile* cmdf = *ciftiMapIter;
        if (cmdf->isEmpty() == false) {
            if (cmdf->isVolumeMappable()) {
                volumeInterfaces.push_back(cmdf);
            }
        }
    }

    const std::vector<int32_t> mapIndices(mapIndicesSet.begin(),
                                          mapIndicesSet.end());
    const int32_t numMapIndices = static_cast<int32_t>(mapIndices.size());
    
    /*
     * In first loop, show values for 'idVolumeFile' (the underlay volume)
     * In second loop, show values for all other volume files
     */
    const VolumeMappableInterface* volumeInterfaceFile = dynamic_cast<VolumeMappableInterface*>(mapFile);
    if (volumeInterfaceFile != NULL) {
            const VolumeFile* volumeFile = dynamic_cast<const VolumeFile*>(volumeInterfaceFile);
            const CiftiMappableDataFile* ciftiFile = dynamic_cast<const CiftiMappableDataFile*>(volumeInterfaceFile);
            CaretAssert((volumeFile != NULL)
                        || (ciftiFile != NULL));
            const CaretMappableDataFile* caretMappableDataFile = dynamic_cast<const CaretMappableDataFile*>(volumeInterfaceFile);
            CaretAssert(caretMappableDataFile != NULL);

            int64_t vfI, vfJ, vfK;
            volumeInterfaceFile->enclosingVoxel(x, y, z,
                                                vfI, vfJ, vfK);

            if (volumeInterfaceFile->indexValid(vfI, vfJ, vfK)) {
                if (volumeFile != NULL) {
                    AString ijkText("IJK ("
                                     + AString::number(vfI)
                                     + ", "
                                     + AString::number(vfJ)
                                     + ", "
                                     + AString::number(vfK)
                                     + ")  ");

                    AString text;
                    
                    for (int32_t k = 0; k < numMapIndices; k++) {
                        CaretAssertVectorIndex(mapIndices, k);
                        const int32_t mapIndex = mapIndices[k];
                        if (k > 0) {
                            text += " ";
                        }
                        if (volumeFile != NULL) {
                            if (volumeFile->getType() == SubvolumeAttributes::LABEL) {
                                const int32_t labelIndex = static_cast<int32_t>(volumeFile->getValue(vfI, vfJ, vfK, mapIndex));
                                const GiftiLabelTable* glt = volumeFile->getMapLabelTable(mapIndex);
                                const GiftiLabel* gl = glt->getLabel(labelIndex);
                                if (gl != NULL) {
                                    text += gl->getName();
                                    text += (" ("
                                             + volumeFile->getMapName(mapIndex)
                                             + ")");
                                }
                                else {
                                    text += ("LABLE_MISSING_FOR_INDEX="
                                             + AString::number(labelIndex));
                                }
                            }
                            else if (volumeFile->getType() == SubvolumeAttributes::RGB) {
                                if (volumeFile->getNumberOfComponents() == 4) {
                                    text += ("RGBA("
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 0))
                                             + ","
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 1))
                                             + ","
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 2))
                                             + ","
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 3))
                                             + ")");
                                }
                                else if (volumeFile->getNumberOfComponents() == 3) {
                                    text += ("RGB("
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 0))
                                             + ","
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 1))
                                             + ","
                                             + AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex, 2))
                                             + ")");
                                }
                            }
                            else {
                                text += AString::number(volumeFile->getValue(vfI, vfJ, vfK, mapIndex));
                            }
                        }
                        else if (ciftiFile != NULL) {

                        }
                    }

                    AString filename;
                    if (dynamic_cast<const VolumeDynamicConnectivityFile*>(volumeFile) != NULL) {
                        filename.append(DataFileTypeEnum::toOverlayTypeName(DataFileTypeEnum::VOLUME_DYNAMIC) + " ");
                    }
                    filename.append(volumeFile->getFileNameNoPath());
                    if (volumeFile->isMappedWithLabelTable()) {
                        labelHtmlTableBuilder.addRow(text,
                                                     //ijkText,
                                                     filename);
                    }
                    else if (volumeFile->isMappedWithPalette()) {
                        scalarHtmlTableBuilder.addRow(text,
                                                      //ijkText,
                                                      filename);
                    }
                }
                else if (ciftiFile != NULL) {
                    if (ciftiFile->isEmpty() == false) {
                        /*
                         * Does file have both label and scalar data
                         */
                        const bool parcelDataFlag(isParcelAndScalarTypeFile(ciftiFile->getDataFileType()));

                        AString textValue;
                        int64_t voxelIJK[3];
                        const QString separator("<br>");
                        if (ciftiFile->getVolumeVoxelIdentificationForMaps(mapIndices,
                                                                           xyz,
                                                                           separator,
                                                                           voxelIJK,
                                                                           textValue)) {
                            AString typeIJKText = (DataFileTypeEnum::toOverlayTypeName(ciftiFile->getDataFileType())
                                                + " "
                                                + "IJK ("
                                                + AString::number(voxelIJK[0])
                                                + ", "
                                                + AString::number(voxelIJK[1])
                                                + ", "
                                                + AString::number(voxelIJK[2])
                                                + ")  ");
                            
                            
                            AString labelText;
                            if (ciftiFile->isMappedWithLabelTable()) {
                                labelText = textValue;
                            }
                            AString scalarText;
                            if (ciftiFile->isMappedWithPalette()) {
                                scalarText = textValue;
                            }
                            if (parcelDataFlag) {
                                /*
                                 * Parcel data has parcel name, separator, scalar value 1, separator, scalar value 2, etc.
                                 * Display parcel name in label table, and scalar
                                 * value in the scalar table
                                 */
                                QStringList parcelAndValue = textValue.split(separator);
                                if (parcelAndValue.size() == static_cast<int>(mapIndices.size() + 1)) {
                                    labelText  = parcelAndValue.at(0);
                                    parcelAndValue.removeAt(0);
                                    scalarText = parcelAndValue.join(separator);
                                }
                            }
                            
                            if ( ! labelText.isEmpty()) {
                                labelHtmlTableBuilder.addRow(labelText,
                                                             //DataFileTypeEnum::toOverlayTypeName(cmdf->getDataFileType()),
                                                             ciftiFile->getFileNameNoPath());
                            }
                            if ( ! scalarText.isEmpty()) {
                                scalarHtmlTableBuilder.addRow(scalarText,
                                                              //DataFileTypeEnum::toOverlayTypeName(cmdf->getDataFileType()),
                                                              ciftiFile->getFileNameNoPath());
                            }
                        }
                    }
                }
            }
        }
}

/**
 * Generate identification text for a surface vertex.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param brain
 *     The brain.
 * @param browserTabContent
 *     Content of the browser tab.
 * @param idSurfaceNode
 *     Information for surface node ID.
 */
void 
IdentificationFormattedTextGenerator::generateSurfaceVertexIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                              const Brain* /*brain*/,
                                                                              const SelectionItemSurfaceNode* idSurfaceNode) const
{
    const Surface* surface = idSurfaceNode->getSurface();
    const int32_t nodeNumber = idSurfaceNode->getNodeNumber();
    
    if ((surface != NULL) 
        && (nodeNumber >= 0)) {
        const float* xyz = surface->getCoordinate(nodeNumber);
        const QString xyzText(AString::number(xyz[0])
                              + ", "
                              + AString::number(xyz[1])
                              + ", "
                              + AString::number(xyz[2]));
        htmlTableBuilder.addHeaderRow(("VERTEX " + QString::number(nodeNumber)),
                                      xyzText,
                                      StructureEnum::toGuiName(surface->getStructure()));
    }
}

/**
 * @return Is the given data file type a parcel and scalar type file?
 * @param dataFileType
 * Type of data file
 */
bool
IdentificationFormattedTextGenerator::isParcelAndScalarTypeFile(const DataFileTypeEnum::Enum dataFileType) const
{
    bool parcelDataFlag      = false;
    
    switch (dataFileType) {
        case DataFileTypeEnum::ANNOTATION:
            break;
        case DataFileTypeEnum::ANNOTATION_TEXT_SUBSTITUTION:
            break;
        case DataFileTypeEnum::BORDER:
            break;
        case DataFileTypeEnum::CONNECTIVITY_DENSE:
            break;
        case DataFileTypeEnum::CONNECTIVITY_DENSE_DYNAMIC:
            break;
        case DataFileTypeEnum::CONNECTIVITY_DENSE_LABEL:
            break;
        case DataFileTypeEnum::CONNECTIVITY_DENSE_PARCEL:
            parcelDataFlag = true;
            break;
        case DataFileTypeEnum::CONNECTIVITY_DENSE_SCALAR:
            break;
        case DataFileTypeEnum::CONNECTIVITY_DENSE_TIME_SERIES:
            break;
        case DataFileTypeEnum::CONNECTIVITY_FIBER_ORIENTATIONS_TEMPORARY:
            break;
        case DataFileTypeEnum::CONNECTIVITY_FIBER_TRAJECTORY_TEMPORARY:
            break;
        case DataFileTypeEnum::CONNECTIVITY_PARCEL:
            parcelDataFlag = true;
            break;
        case DataFileTypeEnum::CONNECTIVITY_PARCEL_DENSE:
            break;
        case DataFileTypeEnum::CONNECTIVITY_PARCEL_LABEL:
            break;
        case DataFileTypeEnum::CONNECTIVITY_PARCEL_SCALAR:
            parcelDataFlag = true;
            break;
        case DataFileTypeEnum::CONNECTIVITY_PARCEL_SERIES:
            parcelDataFlag = true;
            break;
        case DataFileTypeEnum::CONNECTIVITY_SCALAR_DATA_SERIES:
            break;
        case DataFileTypeEnum::FOCI:
            break;
        case DataFileTypeEnum::IMAGE:
            break;
        case DataFileTypeEnum::LABEL:
            break;
        case DataFileTypeEnum::METRIC:
            break;
        case DataFileTypeEnum::METRIC_DYNAMIC:
            break;
        case DataFileTypeEnum::PALETTE:
            break;
        case DataFileTypeEnum::RGBA:
            break;
        case DataFileTypeEnum::SCENE:
            break;
        case DataFileTypeEnum::SPECIFICATION:
            break;
        case DataFileTypeEnum::SURFACE:
            break;
        case DataFileTypeEnum::UNKNOWN:
            CaretAssert(0);
            break;
        case DataFileTypeEnum::VOLUME:
            break;
        case DataFileTypeEnum::VOLUME_DYNAMIC:
            break;
    }

    return parcelDataFlag;
}

/**
 * Generate identification text for a surface vertex.
 * @param labelHtmlTableBuilder
 *     HTML table builder for label identification text.
 * @param scalarHtmlTableBuilder
 *     HTML table builder for scalar identification text.
 * @param mapFile
 *     FIle for generating identification
 * @param mapIndices
 *     Indices of map
 * @param brain
 *     The brain.
 * @param browserTabContent
 *     Content of the browser tab.
 * @param idSurfaceNode
 *     Information for surface node ID.
 */
void
IdentificationFormattedTextGenerator::generateSurfaceDataIdentificationText(HtmlTableBuilder& labelHtmlTableBuilder,
                                                                            HtmlTableBuilder& scalarHtmlTableBuilder,
                                                                            CaretMappableDataFile* mapFile,
                                                                            const std::set<int32_t>& mapIndicesSet,
                                                                            const Brain* /*brain*/,
                                                                            const SelectionItemSurfaceNode* idSurfaceNode) const
{
    const Surface* surface = idSurfaceNode->getSurface();
    const int32_t nodeNumber = idSurfaceNode->getNodeNumber();
    
    if ((surface != NULL)
        && (nodeNumber >= 0)) {
        const BrainStructure* brainStructure = surface->getBrainStructure();
        CaretAssert(brainStructure);
        
        LabelFile* labelFile(NULL);
        MetricFile* metricFile(NULL);
        MetricDynamicConnectivityFile* metricDynConnFile(NULL);
        CiftiMappableDataFile* cmdf = dynamic_cast<CiftiMappableDataFile*>(mapFile);
        if (cmdf == NULL) {
            labelFile = dynamic_cast<LabelFile*>(mapFile);
            if (labelFile == NULL) {
                metricDynConnFile = dynamic_cast<MetricDynamicConnectivityFile*>(mapFile);
                if (metricDynConnFile == NULL) {
                    metricFile = dynamic_cast<MetricFile*>(mapFile);
                }
            }
        }

        std::vector<int32_t> mapIndices(mapIndicesSet.begin(),
                                        mapIndicesSet.end());
        const int numMapIndices = static_cast<int32_t>(mapIndices.size());

        if (cmdf != NULL) {
            AString boldText = (DataFileTypeEnum::toOverlayTypeName(cmdf->getDataFileType())
                                + " "
                                + cmdf->getFileNameNoPath());
            
            /*
             * Does file have both label and scalar data
             */
            const bool parcelDataFlag(isParcelAndScalarTypeFile(cmdf->getDataFileType()));
            AString textValue;
            
            const AString separator("<br>");
            const bool valid = cmdf->getSurfaceNodeIdentificationForMaps(mapIndices,
                                                                         surface->getStructure(),
                                                                         nodeNumber,
                                                                         surface->getNumberOfNodes(),
                                                                         separator,
                                                                         textValue);
            if (valid) {
                AString labelText;
                if (cmdf->isMappedWithLabelTable()) {
                    labelText = textValue;
                }
                AString scalarText;
                if (cmdf->isMappedWithPalette()) {
                    scalarText = textValue;
                }
                if (parcelDataFlag) {
                    /*
                     * Parcel data has parcel name, separator, scalar value 1, separator, scalar value 2, etc.
                     * Display parcel name in label table, and scalar
                     * value in the scalar table
                     */
                    QStringList parcelAndValue = textValue.split(separator);
                    if (parcelAndValue.size() == static_cast<int>(mapIndices.size() + 1)) {
                        labelText  = parcelAndValue.at(0);
                        parcelAndValue.removeAt(0);
                        scalarText = parcelAndValue.join(separator);
                    }
                }
                
                if ( ! labelText.isEmpty()) {
                    labelHtmlTableBuilder.addRow(labelText,
                                                 //DataFileTypeEnum::toOverlayTypeName(cmdf->getDataFileType()),
                                                 cmdf->getFileNameNoPath());
                }
                if ( ! scalarText.isEmpty()) {
                    scalarHtmlTableBuilder.addRow(scalarText,
                                                  //DataFileTypeEnum::toOverlayTypeName(cmdf->getDataFileType()),
                                                  cmdf->getFileNameNoPath());
                }
            }
        }
        
        if (labelFile != NULL) {
            AString text;
            for (int32_t k = 0; k < numMapIndices; k++) {
                CaretAssertVectorIndex(mapIndices, k);
                const int32_t mapIndex = mapIndices[k];
                if (k >= 1) {
                    text.append("<br>");
                }
                AString labelName = labelFile->getLabelName(nodeNumber, mapIndex);
                if (labelName.isEmpty()) {
                    labelName = ("Map-" + AString::number(mapIndex + 1));
                }
                text.append(labelName);
                text.append(" ("
                            + labelFile->getMapName(mapIndex)
                            + ")");
            }
            labelHtmlTableBuilder.addRow(text,
                                         //"LABEL",
                                         labelFile->getFileNameNoPath());
        }
        
        if (metricFile != NULL) {
            AString text;
            for (int32_t k = 0; k < numMapIndices; k++) {
                CaretAssertVectorIndex(mapIndices, k);
                const int32_t mapIndex = mapIndices[k];
                if (k >= 1) {
                    text.append("<br>");
                }
                text.append(AString::number(metricFile->getValue(nodeNumber, mapIndex)));
            }
            scalarHtmlTableBuilder.addRow(text,
                                          //"METRIC",
                                          metricFile->getFileNameNoPath());
        }
        if (metricDynConnFile != NULL) {
            if (metricDynConnFile->isDataValid()) {
                if (metricDynConnFile->isEnabledAsLayer()) {
                    AString boldText = "METRIC DYNAMIC " + metricDynConnFile->getFileNameNoPath();
                    AString text;
                    for (int32_t k = 0; k < numMapIndices; k++) {
                        CaretAssertVectorIndex(mapIndices, k);
                        const int32_t mapIndex = mapIndices[k];
                        if (k >= 1) {
                            text.append("<br>");
                        }
                        text += (" " + AString::number(metricDynConnFile->getValue(nodeNumber, mapIndex)));
                    }
                    scalarHtmlTableBuilder.addRow(text,
                                                  //"METRIC DYNAMIC",
                                                  metricDynConnFile->getFileNameNoPath());
                }
            }
        }
    }
}

/**
 * Find the usage of the file's maps in all overlays.
 *
 * @param caretMappableDataFile
 *    The file whose usage is desired.
 * @param mapIndicesOut
 *    Indices of maps of the file that are used in overlays.
 */
void
IdentificationFormattedTextGenerator::getMapIndicesOfFileUsedInOverlays(const CaretMappableDataFile* caretMappableDataFile,
                                                               std::vector<int32_t>& mapIndicesOut) const
{
    mapIndicesOut.clear();
    
    EventBrowserTabGetAll allTabsEvent;
    EventManager::get()->sendEvent(allTabsEvent.getPointer());
    const std::vector<BrowserTabContent*> allTabs = allTabsEvent.getAllBrowserTabs();
    for (std::vector<BrowserTabContent*>::const_iterator tabIter = allTabs.begin();
         tabIter != allTabs.end();
         tabIter++) {
        BrowserTabContent* tabContent = *tabIter;
        OverlaySet* overlaySet = tabContent->getOverlaySet();
        if (overlaySet != NULL) {
            std::vector<int32_t> mapIndices;
            overlaySet->getSelectedMapIndicesForFile(caretMappableDataFile,
                                                     false,  // true => enabled overlays
                                                     mapIndices);
            mapIndicesOut.insert(mapIndicesOut.end(),
                                 mapIndices.begin(),
                                 mapIndices.end());
        }
    }
    
    /*
     * Sort and remove all duplicates
     */
    if (mapIndicesOut.empty() == false) {
        std::sort(mapIndicesOut.begin(),
                  mapIndicesOut.end());
        std::vector<int32_t>::iterator uniqueIter = std::unique(mapIndicesOut.begin(),
                                                                mapIndicesOut.end());
        mapIndicesOut.resize(std::distance(mapIndicesOut.begin(),
                                        uniqueIter));
    }
}

/**
 * Generate identification text for a data series chart.
 * @param htmlTableBuilder
 *     String builder for identification text.
 * @param idChartDataSeries
 *     Information for chart id.
 */
void
IdentificationFormattedTextGenerator::generateChartDataSeriesIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                       const SelectionItemChartDataSeries* idChartDataSeries) const
{
    if (idChartDataSeries->isValid()) {
        const ChartDataCartesian* chartDataCartesian = idChartDataSeries->getChartDataCartesian();
        
        const ChartDataSource* chartDataSource = chartDataCartesian->getChartDataSource();
        generateChartDataSourceText(htmlTableBuilder,
                                    "DATA SERIES CHART",
                                    chartDataSource);
    }
}

/**
 * Generate identification text for a data series chart.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idChartDataSeries
 *     Information for chart id.
 */
void
IdentificationFormattedTextGenerator::generateChartFrequencySeriesIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                       const SelectionItemChartFrequencySeries* idChartFrequencySeries) const
{
    if (idChartFrequencySeries->isValid()) {
        const ChartDataCartesian* chartDataCartesian = idChartFrequencySeries->getChartDataCartesian();
        
        const ChartDataSource* chartDataSource = chartDataCartesian->getChartDataSource();
        generateChartDataSourceText(htmlTableBuilder,
                                    "FREQUENCY SERIES CHART",
                                    chartDataSource);
    }
}

/**
 * Generate identification text for a matrix chart.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idChartMatrix
 *     Information for matrix chart id.
 */
void
IdentificationFormattedTextGenerator::generateChartMatrixIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                   const SelectionItemChartMatrix* idChartMatrix) const
{
    if (idChartMatrix->isValid()) {
        const ChartableMatrixInterface* chartMatrixInterface = idChartMatrix->getChartableMatrixInterface();
        const CaretMappableDataFile* caretMappableDataFile = chartMatrixInterface->getMatrixChartCaretMappableDataFile();
        
        const int32_t rowIndex = idChartMatrix->getMatrixRowIndex();
        const int32_t columnIndex = idChartMatrix->getMatrixColumnIndex();
        AString rowName;
        AString columnName;
        AString cellValue;
        const bool validData = chartMatrixInterface->getMatrixCellAttributes(rowIndex,
                                                                             columnIndex,
                                                                             cellValue,
                                                                             rowName,
                                                                             columnName);
        
        AString boldText("MATRIX CHART");
        
        if (validData) {
            htmlTableBuilder.addRow(cellValue,
                                    (("Row: " + AString::number(rowIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI()))
                                     + "<br>"
                                     + ("Column: " + AString::number(columnIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI()))),
                                    caretMappableDataFile->getFileNameNoPath());
        }
        else {
            htmlTableBuilder.addRow(boldText,
                                    "",
                                    caretMappableDataFile->getFileNameNoPath());
        }
    }
}

/**
 * Generate identification text for a chart two histogram.
 *
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idText
 *     Identification string builder
 * @param idChartTwoHistogram
 *     Information for selected chart two histogram.
 * @param mapFile
 *     FIle for generating identification
 * @param mapIndices
 *     Indices of map
 * @param toolTipFlag
 *     If true, create tooltip text
 */
void
IdentificationFormattedTextGenerator::generateChartTwoHistogramIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                                  IdentificationStringBuilder& idText,
                                                                                  const SelectionItemChartTwoHistogram* idChartTwoHistogram,
                                                                                  CaretMappableDataFile* mapFile,
                                                                                  const std::set<int32_t>& mapIndices,
                                                                                  const bool toolTipFlag) const
{
    const int32_t mapIndex    = idChartTwoHistogram->getMapIndex();
    const int32_t bucketIndex = idChartTwoHistogram->getBucketIndex();
    const bool    allMapsFlag = idChartTwoHistogram->isAllMapsSelected();
    
    if (idChartTwoHistogram->isValid()) {
        ChartableTwoFileHistogramChart* fileHistogramChart = idChartTwoHistogram->getFileHistogramChart();
        CaretAssert(fileHistogramChart);
        CaretMappableDataFile* chartMapFile = fileHistogramChart->getCaretMappableDataFile();
        CaretAssert(chartMapFile);
        
        {
            ChartableTwoFileHistogramChart* chartingDelegate = chartMapFile->getChartingDelegate()->getHistogramCharting();
            CaretAssert(chartingDelegate);
            const Histogram* histogram = chartingDelegate->getHistogramForChartDrawing(mapIndex,
                                                                                       allMapsFlag);
            CaretAssert(histogram);
            
            float bucketValue = 0.0;
            float bucketHeight = 0.0;
            QString columnOne;
            if (histogram->getHistogramDisplayBucketDataValueAndHeight(bucketIndex, bucketValue, bucketHeight)) {
                if (toolTipFlag) {
                    AString boldText("Histogram");
                    idText.addLine(false,
                                   boldText,
                                   chartMapFile->getFileNameNoPath());
                    
                    idText.addLine(true,
                                   "Bucket Index",
                                   (AString::number(bucketIndex)));
                    
                    idText.addLine(true,
                                   "Data Value at Bucket",
                                   (AString::number(bucketValue)));
                    
                    const int64_t bucketHeightInteger = static_cast<int64_t>(bucketHeight);
                    idText.addLine(true,
                                   "Bucket Count",
                                   (AString::number(bucketHeightInteger)));
                }
                else {
                    if (mapFile == chartMapFile) {
                        if ((mapIndices.find(mapIndex) != mapIndices.end())
                            || allMapsFlag) {
                            columnOne.append("Bucket: " + (AString::number(bucketIndex)));
                            columnOne.append("<br>Data Value: " +  (AString::number(bucketValue)));
                            
                            const int64_t bucketHeightInteger = static_cast<int64_t>(bucketHeight);
                            columnOne.append("<br>Bucket Count: " + (AString::number(bucketHeightInteger)));
                            
                            htmlTableBuilder.addRow(columnOne,
                                                    "Histogram",
                                                    chartMapFile->getFileNameNoPath());
                        }
                    }
                }
            }
        }
    }
}

/**
 * Generate identification text for a chart two line-series.
 *
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idText
 *     Identification string builder
 * @param idChartTwoLineSeries
 *     Information for selected chart two line-series.
 * @param mapFile
 *     FIle for generating identification
 * @param mapIndices
 *     Indices of map
 * @param toolTipFlag
 *     If true, generate tooltip
 */
void
IdentificationFormattedTextGenerator::generateChartTwoLineSeriesIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                                   IdentificationStringBuilder& idText,
                                                                                   const SelectionItemChartTwoLineSeries* idChartTwoLineSeries,
                                                                                   CaretMappableDataFile* mapFile,
                                                                                   const std::set<int32_t>& /*mapIndices*/,
                                                                                   const bool toolTipFlag) const
{
    if (idChartTwoLineSeries->isValid()) {
        const ChartableTwoFileLineSeriesChart* fileLineSeriesChart = idChartTwoLineSeries->getFileLineSeriesChart();
        CaretAssert(fileLineSeriesChart);
        const CaretMappableDataFile* chartMapFile = fileLineSeriesChart->getCaretMappableDataFile();
        CaretAssert(chartMapFile);
        const ChartTwoDataCartesian* cartesianData = idChartTwoLineSeries->getChartTwoCartesianData();
        CaretAssert(cartesianData);
        const MapFileDataSelector* mapFileDataSelector = cartesianData->getMapFileDataSelector();
        CaretAssert(mapFileDataSelector);
        
        if ( ! toolTipFlag) {
            if (chartMapFile != mapFile) {
                return;
            }
        }
        int32_t primitiveIndex = idChartTwoLineSeries->getLineSegmentIndex();
        
        AString boldText("Line Chart");
        
        cartesianData->getGraphicsPrimitive();
        const GraphicsPrimitive* primitive = cartesianData->getGraphicsPrimitive();
        CaretAssert(primitive);
        
        if (primitiveIndex >= 1) {
            float xyz1[3];
            primitive->getVertexFloatXYZ(primitiveIndex - 1,
                                         xyz1);
            float xyz2[3];
            primitive->getVertexFloatXYZ(primitiveIndex,
                                         xyz2);
            if (toolTipFlag) {
                idText.addLine(true,
                               "XY Start",
                               AString::fromNumbers(xyz1, 2, ", "));
                idText.addLine(true,
                               "XY End ",
                               AString::fromNumbers(xyz2, 2, ", "));
            }
            else {
                htmlTableBuilder.addRow(("XY Start:" + AString::fromNumbers(xyz1, 2, ", "))
                                        + ("XY End:" + AString::fromNumbers(xyz2, 2, ", ")),
                                        boldText,
                                        chartMapFile->getFileNameNoPath());
            }
        }
        else {
            float xyz[3];
            primitive->getVertexFloatXYZ(primitiveIndex,
                                         xyz);
            if (toolTipFlag) {
                idText.addLine(true,
                               "XY",
                               AString::fromNumbers(xyz, 2, ", "));
            }
            else {
                htmlTableBuilder.addRow(("XY:" + AString::fromNumbers(xyz, 2, ", ")),
                                        boldText,
                                        chartMapFile->getFileNameNoPath());
            }
        }
        
        generateMapFileSelectorText(htmlTableBuilder,
                                    mapFileDataSelector);
    }
}

/**
 * Generate identification text for a chart two matrix.
 *
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idText
 *     Identification string builder
 * @param idChartTwoMatrix
 *     Information for selected chart two matrix.
 * @param mapFile
 *     FIle for generating identification
 * @param mapIndices
 *     Indices of map
 * @param toolTipFlag
 *     If true, generate tool tip
 */
void
IdentificationFormattedTextGenerator::generateChartTwoMatrixIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                               IdentificationStringBuilder& idText,
                                                                               const SelectionItemChartTwoMatrix* idChartTwoMatrix,
                                                                               CaretMappableDataFile* mapFile,
                                                                               const std::set<int32_t>& /*mapIndices*/,
                                                                               const bool toolTipFlag) const
{
    if (idChartTwoMatrix->isValid()) {
        const ChartableTwoFileMatrixChart* matrixChart = idChartTwoMatrix->getFileMatrixChart();
        CaretAssert(matrixChart);
        
        const int32_t rowIndex = idChartTwoMatrix->getRowIndex();
        const int32_t colIndex = idChartTwoMatrix->getColumnIndex();
        
        const CaretMappableDataFile* chartMapFile = matrixChart->getCaretMappableDataFile();
        CaretAssert(chartMapFile);
       
        if ( ! toolTipFlag) {
            if (chartMapFile != mapFile) {
                return;
            }
        }
        
        AString boldText("MATRIX ");
        QString rowText;
        QString colText;
        if ((rowIndex >= 0)
            && (matrixChart->hasRowSelection())) {
            const AString rowName = matrixChart->getRowName(rowIndex);
            rowText.append("Row "
                           + AString::number(rowIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI())
                           + ": "
                           + rowName);
        }
        if ((colIndex >= 0)
            && (matrixChart->hasColumnSelection())) {
            const AString colName = matrixChart->getColumnName(colIndex);
            colText.append("Column "
                           + AString::number(colIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI())
                           + ": "
                           + colName);
        }
        
        if (( ! colText.isEmpty())
            || ( ! rowText.isEmpty())) {
            if (toolTipFlag) {
                if ( ! colText.isEmpty()) {
                    colText.append(" ");
                }
                idText.addLine(true,
                               (colText + rowText),
                               chartMapFile->getFileNameNoPath());
            }
            else {
                htmlTableBuilder.addRow(rowText,
                                        colText,
                                        chartMapFile->getFileNameNoPath());
            }
        }
    }
}

/**
 * Generate identification text for a CIFTI Connectivity Matrix Row/Column
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idCiftiConnMatrix
 *     Information for CIFTI Connectivity Matrix Row/Column.
 */
void
IdentificationFormattedTextGenerator::generateCiftiConnectivityMatrixIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                               const SelectionItemCiftiConnectivityMatrixRowColumn* idCiftiConnMatrix) const
{
    if (idCiftiConnMatrix->isValid()) {
        const CiftiMappableConnectivityMatrixDataFile* connMatrixFile = idCiftiConnMatrix->getCiftiConnectivityMatrixFile();
        const int32_t rowIndex = idCiftiConnMatrix->getMatrixRowIndex();
        const int32_t colIndex = idCiftiConnMatrix->getMatrixColumnIndex();
        
        AString colTwoText("MATRIX ROW/COLUMN");
        AString colThreeText(connMatrixFile->getFileNameNoPath());
        
        AString colOneText;
        AString rowName = " ";
        AString colName = " ";
        bool validData = true;
        if (validData) {
            if (rowIndex >= 0) {
                colOneText = ("Row " + AString::number(rowIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI())
                              + "<br>"
                              + rowName);
            }
            if (colIndex >= 0) {
                colOneText = ("Column " + AString::number(colIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI())
                              + "<br>"
                              + colName);
            }
        }
        
        htmlTableBuilder.addRow(colOneText,
                                colTwoText,
                                colThreeText);
    }
}

/**
 * Generate identification text for chart data source.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param typeOfChartText
 *     Text describing the type of chart.
 * @param chartDataSource
 *     Source of chart data.
 */
void
IdentificationFormattedTextGenerator::generateChartDataSourceText(HtmlTableBuilder& htmlTableBuilder,
                                                                  const AString& typeOfChartText,
                                                                  const ChartDataSource* chartDataSource) const
{
    AString chartFileName = chartDataSource->getChartableFileName();
    if (! chartFileName.isEmpty()) {
        chartFileName = FileInformation(chartFileName).getFileName();
    }
    
    AString columnOne;
    switch (chartDataSource->getDataSourceMode()) {
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_INVALID:
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_FILE_ROW:
        {
            AString fileName;
            int32_t rowIndex;
            chartDataSource->getFileRow(fileName,
                                        rowIndex);
            chartFileName = fileName;
            columnOne = ("Row " + AString::number(rowIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI()));
        }
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDEX:
        {
            AString structureName;
            int32_t numberOfNodes;
            int32_t nodeIndex;
            chartDataSource->getSurfaceNode(structureName,
                                            numberOfNodes,
                                            nodeIndex);
            columnOne = ("Vertex Index " + AString::number(nodeIndex));
            columnOne.append("<br>Structure " + structureName);
        }
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDICES_AVERAGE:
        {
            AString structureName;
            int32_t numberOfNodes;
            std::vector<int32_t> nodeIndices;
            chartDataSource->getSurfaceNodeAverage(structureName, numberOfNodes, nodeIndices);
            columnOne = ("Vertex Avg Count " + AString::number(nodeIndices.size()));
            columnOne.append("<br>Structure " + structureName);
        }
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_VOXEL_IJK:
        {
            float voxelXYZ[3];
            chartDataSource->getVolumeVoxel(voxelXYZ);
            columnOne = ("Voxel XYZ " + AString::fromNumbers(voxelXYZ, 3, ","));
        }
            break;
    }
    
    if (! columnOne.isEmpty()) {
        htmlTableBuilder.addRow(columnOne,
                                typeOfChartText,
                                chartFileName);
    }
}

/**
 * Generate text for a map file data selector.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param mapFileDataSelector
 *     The map file data selector.
 */
void
IdentificationFormattedTextGenerator::generateMapFileSelectorText(HtmlTableBuilder& htmlTableBuilder,
                                                                  const MapFileDataSelector* mapFileDataSelector) const
{
    
    switch (mapFileDataSelector->getDataSelectionType()) {
        case MapFileDataSelector::DataSelectionType::INVALID:
            break;
        case MapFileDataSelector::DataSelectionType::COLUMN_DATA:
        {
            CaretMappableDataFile* mapFile = NULL;
            AString mapFileName;
            int32_t columnIndex = -1;
            mapFileDataSelector->getColumnIndex(mapFile,
                                                mapFileName,
                                                columnIndex);
            htmlTableBuilder.addRow("Column: " + AString::number(columnIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI()),
                                    "",
                                    mapFileName);
        }
            break;
        case MapFileDataSelector::DataSelectionType::ROW_DATA:
        {
            CaretMappableDataFile* mapFile = NULL;
            AString mapFileName;
            int32_t rowIndex = -1;
            mapFileDataSelector->getRowIndex(mapFile,
                                             mapFileName,
                                             rowIndex);
            htmlTableBuilder.addRow("Row: " + AString::number(rowIndex + CiftiMappableDataFile::getCiftiFileRowColumnIndexBaseForGUI()),
                                    "",
                                    mapFileName);
        }
            break;
        case MapFileDataSelector::DataSelectionType::SURFACE_VERTEX:
        {
            StructureEnum::Enum structure = StructureEnum::INVALID;
            int32_t numberOfVertices = 0;
            int32_t vertexIndex = -1;
            mapFileDataSelector->getSurfaceVertex(structure,
                                                  numberOfVertices,
                                                  vertexIndex);
            
            if ((structure != StructureEnum::INVALID)
                && (vertexIndex >= 0)) {
                htmlTableBuilder.addRow("Vertex: " + AString::number(vertexIndex),
                                        "Structure: " + StructureEnum::toGuiName(structure),
                                        "");
            }
        }
            break;
        case MapFileDataSelector::DataSelectionType::SURFACE_VERTICES_AVERAGE:
        {
            StructureEnum::Enum structure = StructureEnum::INVALID;
            int32_t numberOfVertices = 0;
            std::vector<int32_t> vertexIndices;
            mapFileDataSelector->getSurfaceVertexAverage(structure,
                                                         numberOfVertices,
                                                         vertexIndices);
            
            const int32_t averageCount = static_cast<int32_t>(vertexIndices.size());
            if ((structure != StructureEnum::INVALID)
                && (averageCount > 0)) {
                htmlTableBuilder.addRow("Vertex Avg Count: " + AString::number(averageCount),
                                        "Structure: " + StructureEnum::toGuiName(structure),
                                        "");
            }
        }
            break;
        case MapFileDataSelector::DataSelectionType::VOLUME_XYZ:
        {
            float voxelXYZ[3];
            mapFileDataSelector->getVolumeVoxelXYZ(voxelXYZ);
            htmlTableBuilder.addRow("Voxel XYZ: " + AString::fromNumbers(voxelXYZ, 3, ","));
        }
            break;
    }
}


/**
 * Generate identification text for a time series chart.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idChartTimeSeries
 *     Information for chart id.
 */
void
IdentificationFormattedTextGenerator::generateChartTimeSeriesIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                                       const SelectionItemChartTimeSeries* idChartTimeSeries) const
{
    if (idChartTimeSeries->isValid()) {
        const ChartDataCartesian* chartDataCartesian = idChartTimeSeries->getChartDataCartesian();
        
        const ChartDataSource* chartDataSource = chartDataCartesian->getChartDataSource();
        generateChartDataSourceText(htmlTableBuilder,
                                    "TIME SERIES CHART",
                                    chartDataSource);
    }
}

/**
 * Generate identification text for a surface border identification.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idText
 *     Identification string builder
 * @param idSurfaceBorder
 *     Information for surface border ID.
 * @param toolTipFlag
 *     True if this is for tooltip.
 */
void
IdentificationFormattedTextGenerator::generateSurfaceBorderIdentifcationText(HtmlTableBuilder& htmlTableBuilder,
                                                                             IdentificationStringBuilder& idText,
                                                                    const SelectionItemBorderSurface* idSurfaceBorder,
                                                                    const bool toolTipFlag) const
{
    if (idSurfaceBorder->isValid()) {
        const Border* border = idSurfaceBorder->getBorder();
        const SurfaceProjectedItem* spi = border->getPoint(idSurfaceBorder->getBorderPointIndex());
        float xyz[3];
        spi->getProjectedPosition(*idSurfaceBorder->getSurface(), xyz, false);

        if (toolTipFlag) {
            bool indentFlag = false;
            idText.addLine(indentFlag,
                           "Border",
                           border->getName());
            indentFlag = true;
            idText.addLine(indentFlag,
                           "XYZ",
                           AString::fromNumbers(xyz, 3, ","));
        }
        else {
            const AString numberIndexText = ("("
                                             + AString::number(idSurfaceBorder->getBorderIndex())
                                             + ","
                                             + AString::number(idSurfaceBorder->getBorderPointIndex()));

            htmlTableBuilder.addRow(AString::fromNumbers(xyz, 3, ","),
                                          ("BORDER " + numberIndexText),
                                          ("Name: " + border->getName()
                                           + "<br>Class: " + border->getClassName()));
            
        }
    }
}

/**
 * Generate identification text for a surface focus identification.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idSurfaceFocus
 *     Information for surface focus ID.
 * @param toolTipFlag
 *     True if this is for tooltip.
 */
void
IdentificationFormattedTextGenerator::generateSurfaceFocusIdentifcationText(HtmlTableBuilder& htmlTableBuilder,
                                                                  const SelectionItemFocusSurface* idSurfaceFocus,
                                                                  const bool toolTipFlag) const
{
    if (idSurfaceFocus->isValid()) {
        const Focus* focus = idSurfaceFocus->getFocus();
        const int32_t projectionIndex = idSurfaceFocus->getFocusProjectionIndex();
        
        IdentificationStringBuilder idText;
        generateFocusIdentifcationText(htmlTableBuilder,
                                       idText,
                                       focus,
                                       idSurfaceFocus->getFocusIndex(),
                                       projectionIndex,
                                       idSurfaceFocus->getSurface(),
                                       toolTipFlag);
    }
}

/**
 * Generate identification text for a surface focus identification.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idText
 *     String builder for identification text.
 * @param focus
 *     The focus.
 * @param focusIndex
 *     Index of focus
 * @param projectionIndex
 *     Index of projection
 * @param surface
 *     Surface for focus (may be NULL)
 * @param toolTipFlag
 *     True if this is for tooltip.
 */
void
IdentificationFormattedTextGenerator::generateFocusIdentifcationText(HtmlTableBuilder& htmlTableBuilder,
                                                                     IdentificationStringBuilder& idText,
                                                                     const Focus* focus,
                                                                     const int32_t focusIndex,
                                                                     const int32_t projectionIndex,
                                                                     const Surface* surface,
                                                                     const bool toolTipFlag) const
{
    const SurfaceProjectedItem* spi = focus->getProjection(projectionIndex);
    float xyzStereo[3];
    spi->getStereotaxicXYZ(xyzStereo);
    const AString stereoXYZText((spi->isStereotaxicXYZValid()
                                 ? AString::fromNumbers(xyzStereo, 3, ",")
                                 : "Invalid"));
    if (toolTipFlag) {
        bool indentFlag = false;
        idText.addLine(indentFlag,
                       "Focus",
                       focus->getName());
        indentFlag = true;
        idText.addLine(indentFlag,
                       "XYZ",
                       stereoXYZText);
    }
    else {
        htmlTableBuilder.addRow(stereoXYZText,
                                      ("FOCUS " + AString::number(focusIndex)),
                                      ("Name: " + focus->getName() + "<br>Class: " + focus->getClassName()));
        float xyzProj[3];
        if (surface != NULL) {
            spi->getProjectedPosition(*surface, xyzProj, false);
        }
        
        
        bool projValid = false;
        AString xyzProjName = "XYZ (Projected)";
        if (spi->getBarycentricProjection()->isValid()) {
            xyzProjName = "(Projected to Triangle)";
            projValid = true;
        }
        else if (spi->getVanEssenProjection()->isValid()) {
            xyzProjName = "(Projected to Edge)";
            projValid = true;
        }
        if (projValid) {
            idText.addLine(true,
                           xyzProjName,
                           xyzProj,
                           3,
                           true);
            htmlTableBuilder.addRow((AString::fromNumbers(xyzProj, 3, ", ") + xyzProjName),
                                    StructureEnum::toGuiName(spi->getStructure()),
                                    "");
        }
        else {
            htmlTableBuilder.addRow("Projection Invalid",
                                    StructureEnum::toGuiName(spi->getStructure()),
                                    "");
        }
        //
        //            const int32_t numberOfProjections = focus->getNumberOfProjections();
        //            for (int32_t i = 0; i < numberOfProjections; i++) {
        //                if (i != projectionIndex) {
        //                    const SurfaceProjectedItem* proj = focus->getProjection(i);
        //                    AString projTypeName = "";
        //                    if (proj->getBarycentricProjection()->isValid()) {
        //                        projTypeName = "Triangle";
        //
        //                    }
        //                    else if (proj->getVanEssenProjection()->isValid()) {
        //                        projTypeName = "Edge";
        //                    }
        //                    if (projTypeName.isEmpty() == false) {
        //                        const AString txt = (StructureEnum::toGuiName(proj->getStructure())
        //                                             + " ("
        //                                             + projTypeName
        //                                             + ")");
        //
        //                        idText.addLine(true,
        //                                       "Ambiguous Projection",
        //                                       txt);
        //                    }
        //                }
        //            }
        //
        //            idText.addLine(true,
        //                           "Area",
        //                           focus->getArea());
        //
        //            idText.addLine(true,
        //                           "Class Name",
        //                           focus->getClassName());
        //
        //            idText.addLine(true,
        //                           "Comment",
        //                           focus->getComment());
        //
        //            idText.addLine(true,
        //                           "Extent",
        //                           focus->getExtent(),
        //                           true);
        //
        //            idText.addLine(true,
        //                           "Geography",
        //                           focus->getGeography());
        //
        //            idText.addLine(true,
        //                           "Region of Interest",
        //                           focus->getRegionOfInterest());
        //
        //            idText.addLine(true,
        //                           "Statistic",
        //                           focus->getStatistic());
    }
}

/**
 * Generate identification text for a volume focus identification.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idVolumeFocus
 *     Information for surface focus ID.
 */
void
IdentificationFormattedTextGenerator::generateVolumeFocusIdentifcationText(HtmlTableBuilder& htmlTableBuilder,
                                                                  const SelectionItemFocusVolume* idVolumeFocus) const
{
    if (idVolumeFocus->isValid()) {
        const Focus* focus = idVolumeFocus->getFocus();
        const SurfaceProjectedItem* spi = focus->getProjection(idVolumeFocus->getFocusProjectionIndex());
        float xyzVolume[3];
        spi->getVolumeXYZ(xyzVolume);
        float xyzStereo[3];
        spi->getStereotaxicXYZ(xyzStereo);
        
        IdentificationStringBuilder idText;
        generateFocusIdentifcationText(htmlTableBuilder,
                                       idText,
                                       focus,
                                       idVolumeFocus->getFocusIndex(),
                                       idVolumeFocus->getFocusProjectionIndex(),
                                       NULL,
                                       false);
    }
}

/**
 * Generate identification text for image identification.
 * @param htmlTableBuilder
 *     HTML table builder for identification text.
 * @param idImage
 *     Information for image ID.
 */
void
IdentificationFormattedTextGenerator::generateImageIdentificationText(HtmlTableBuilder& htmlTableBuilder,
                                                             const SelectionItemImage* idImage) const
{
    if (idImage->isValid()) {
        uint8_t pixelRGBA[4] = { 0, 0, 0, 0 };
        idImage->getPixelRGBA(pixelRGBA);
        htmlTableBuilder.addRow((" RGBA (" + AString::fromNumbers(pixelRGBA, 4, ",") + ")"),
                                ("Pixel IJ ("
                                 + AString::number(idImage->getPixelI())
                                 + ","
                                 + AString::number(idImage->getPixelJ())
                                 + ")"),
                                idImage->getImageFile()->getFileNameNoPath());
                                
    }
}

/**
 * Get text for the tooltip for a selected node.
 *
 * @param brain
 *     The Brain.
 * @param browserTab
 *     Browser tab in which tooltip is displayed
 * @param selectionManager
 *     The selection manager.
 * @param dataToolTipsManager
 *     The data tooltips manager
 * @param idText
 *     String builder for identification text.
 */
void
IdentificationFormattedTextGenerator::generateSurfaceToolTip(const Brain* brain,
                                                          const BrowserTabContent* browserTab,
                                                          const SelectionManager* selectionManager,
                                                          const DataToolTipsManager* dataToolTipsManager,
                                                          IdentificationStringBuilder& idText) const
{
    std::unique_ptr<HtmlTableBuilder> htmlTableBuilder = createHtmlTableBuilder(3);
    
    const SelectionItemSurfaceNode* nodeSelection = selectionManager->getSurfaceNodeIdentification();
    CaretAssert(nodeSelection);
    if (nodeSelection->isValid()) {
        const Surface* surface = nodeSelection->getSurface();
        CaretAssert(surface);
        int32_t surfaceNumberOfNodes = surface->getNumberOfNodes();
        int32_t surfaceNodeIndex = nodeSelection->getNodeNumber();
        StructureEnum::Enum surfaceStructure = surface->getStructure();
        
        bool indentFlag = false;
        if ((surfaceStructure != StructureEnum::INVALID)
            && (surfaceNumberOfNodes > 0)
            && (surfaceNodeIndex >= 0)) {
            
            bool addVertexFlag(false);
            bool showSurfaceFlag = dataToolTipsManager->isShowSurfaceViewed();
            if (dataToolTipsManager->isShowSurfacePrimaryAnatomical()) {
                const Surface* anatSurface = brain->getPrimaryAnatomicalSurfaceForStructure(surfaceStructure);
                if (anatSurface != NULL) {
                    if (anatSurface->getNumberOfNodes() == surfaceNumberOfNodes) {
                        float xyz[3];
                        anatSurface->getCoordinate(surfaceNodeIndex,
                                                   xyz);
                        idText.addLine(indentFlag,
                                       "Vertex",
                                       AString::number(surfaceNodeIndex));
                        indentFlag = true;
                        addVertexFlag = false;
                        
                        idText.addLine(indentFlag,
                                       "Anatomy Surface",
                                       AString::fromNumbers(xyz, 3, ", ", 'f', 2));
                        if (surface == anatSurface) {
                            showSurfaceFlag = false;
                        }
                    }
                }
            }
            
            if (showSurfaceFlag) {
                float xyz[3];
                surface->getCoordinate(surfaceNodeIndex,
                                       xyz);
                if (addVertexFlag) {
                    idText.addLine(indentFlag,
                                   "Vertex",
                                   AString::number(surfaceNodeIndex));
                    indentFlag = true;
                }
                
                idText.addLine(indentFlag,
                               (SurfaceTypeEnum::toGuiName(surface->getSurfaceType())
                                + " Surface"),
                               AString::fromNumbers(xyz, 3, ", "));
            }
            
            if (dataToolTipsManager->isShowTopEnabledLayer()) {
                OverlaySet* overlaySet = const_cast<OverlaySet*>(browserTab->getOverlaySet());
                CaretAssert(overlaySet);
                Overlay* overlay = getTopEnabledOverlay(overlaySet);
                if (overlay != NULL) {
                    CaretMappableDataFile* mapFile(NULL);
                    int32_t mapIndex(-1);
                    overlay->getSelectionData(mapFile,
                                              mapIndex);
                    if ((mapFile != NULL)
                        && (mapIndex >= 0)) {
                        std::vector<int32_t> mapIndices { mapIndex };
                        AString textValue;
                        mapFile->getSurfaceNodeIdentificationForMaps(mapIndices,
                                                                     surfaceStructure,
                                                                     surfaceNodeIndex,
                                                                     surfaceNumberOfNodes,
                                                                     " ",
                                                                     textValue);
                        if ( ! textValue.isEmpty()) {
                            idText.addLine(indentFlag,
                                           "Top Enabled Layer",
                                           textValue);
                        }
                    }
                }
            }
        }
    }
    
    if (dataToolTipsManager->isShowBorder()) {
        const SelectionItemBorderSurface* borderSelection = selectionManager->getSurfaceBorderIdentification();
        CaretAssert(borderSelection);
        if (borderSelection->isValid()) {
            generateSurfaceBorderIdentifcationText(*htmlTableBuilder,
                                                   idText,
                                                   borderSelection,
                                                   true);
        }
    }
    
    if (dataToolTipsManager->isShowFocus()) {
        const SelectionItemFocusSurface* focusSelection = selectionManager->getSurfaceFocusIdentification();
        CaretAssert(focusSelection);
        if (focusSelection->isValid()) {
            const FociFile* fociFile = focusSelection->getFociFile();
            const int32_t focusIndex = focusSelection->getFocusIndex();
            if ((fociFile != NULL)
                && (focusIndex >= 0)) {
                const Focus* focus = fociFile->getFocus(focusIndex);
                if (focus != NULL) {
                    generateFocusIdentifcationText(*htmlTableBuilder,
                                                   idText,
                                                   focus,
                                                   focusIndex,
                                                   focusSelection->getFocusProjectionIndex(),
                                                   focusSelection->getSurface(),
                                                   true);
                }
            }
        }
    }
}

/**
 * Get text for the tooltip for a selected node.
 *
 * @param browserTab
 *     Browser tab in which tooltip is displayed
 * @param selectionManager
 *     The selection manager.
 * @param dataToolTipsManager
 *     The data tooltips manager
 * @param idText
 *     String builder for identification text.
 */
void
IdentificationFormattedTextGenerator::generateVolumeToolTip(const BrowserTabContent* browserTab,
                                                   const SelectionManager* selectionManager,
                                                   const DataToolTipsManager* dataToolTipsManager,
                                                   IdentificationStringBuilder& idText) const
{
    const SelectionItemVoxel* voxelSelection = selectionManager->getVoxelIdentification();
    
    OverlaySet* overlaySet = const_cast<OverlaySet*>(browserTab->getOverlaySet());
    CaretAssert(overlaySet);
    
    double selectionXYZ[3];
    voxelSelection->getModelXYZ(selectionXYZ);
    float xyz[3] {
        static_cast<float>(selectionXYZ[0]),
        static_cast<float>(selectionXYZ[1]),
        static_cast<float>(selectionXYZ[2])
    };
    
    bool indentFlag = false;
    if (dataToolTipsManager->isShowVolumeUnderlay()) {
        Overlay* volumeUnderlay = overlaySet->getUnderlayContainingVolume();
        if (volumeUnderlay != NULL) {
            CaretMappableDataFile* mapFile = NULL;
            int32_t mapIndex(-1);
            volumeUnderlay->getSelectionData(mapFile,
                                             mapIndex);
            
            VolumeMappableInterface* underlayVolumeInterface = NULL;
            if (mapFile != NULL) {
                underlayVolumeInterface = dynamic_cast<VolumeMappableInterface*>(mapFile);
                CaretAssert(underlayVolumeInterface == overlaySet->getUnderlayVolume());
            }
            
            if (underlayVolumeInterface != NULL) {
                /*
                 * Update IJK and XYZ since selection XYZ may be
                 * a different volume file.
                 */
                int64_t selectionIJK[3];
                voxelSelection->getVoxelIJK(selectionIJK);
                int64_t ijk[3] { selectionIJK[0], selectionIJK[1], selectionIJK[2] };
                
                
                bool validFlag(false);
                const float value = underlayVolumeInterface->getVoxelValue(xyz[0], xyz[1], xyz[2],
                                                                           &validFlag,
                                                                           mapIndex);
                if (validFlag) {
                    underlayVolumeInterface->enclosingVoxel(xyz[0], xyz[1], xyz[2],
                                                            ijk[0], ijk[1], ijk[2]);
                    underlayVolumeInterface->indexToSpace(ijk, xyz);
                    idText.addLine(indentFlag,
                                   "Underlay Value",
                                   AString::number(value, 'f'));
                    indentFlag = true;
                    idText.addLine(indentFlag,
                                   "IJK: ",
                                   AString::fromNumbers(ijk, 3, ", "));
                    idText.addLine(indentFlag,
                                   "XYZ",
                                   AString::fromNumbers(xyz, 3, ", ", 'f', 1));
                }
            }
        }
    }
    
    if (dataToolTipsManager->isShowTopEnabledLayer()) {
        Overlay* overlay = getTopEnabledOverlay(overlaySet);
        if (overlay != NULL) {
            CaretMappableDataFile* mapFile(NULL);
            int32_t mapIndex(-1);
            overlay->getSelectionData(mapFile,
                                      mapIndex);
            if ((mapFile != NULL)
                && (mapIndex >= 0)) {
                std::vector<int32_t> mapIndices { mapIndex };
                AString textValue;
                int64_t ijk[3];
                mapFile->getVolumeVoxelIdentificationForMaps(mapIndices,
                                                             xyz,
                                                             " ",
                                                             ijk,
                                                             textValue);
                if ( ! textValue.isEmpty()) {
                    idText.addLine(indentFlag,
                                   ("Top Enabled Layer: "
                                    + textValue));
                }
            }
        }
    }
}

/**
 * @return Get the top-most enabled overlay.  NULL if no overlays enabled
 *
 * @param overlaySet
 *     Overlay set for overlay.
 */
Overlay*
IdentificationFormattedTextGenerator::getTopEnabledOverlay(OverlaySet* overlaySet) const
{
    CaretAssert(overlaySet);
    const int32_t numberOfOverlays = overlaySet->getNumberOfDisplayedOverlays();
    for (int32_t i = 0; i < numberOfOverlays; i++) {
        Overlay* overlay = overlaySet->getOverlay(i);
        CaretAssert(overlay);
        if (overlay->isEnabled()) {
            return overlay;
        }
    }
    return NULL;
}

/**
 * Get text for the tooltip for a selected node.
 *
 * @param selectionManager
 *     The selection manager.
 * @param dataToolTipsManager
 *     The data tooltips manager
 * @param idText
 *     String builder for identification text.
 */
void
IdentificationFormattedTextGenerator::generateChartToolTip(const SelectionManager* selectionManager,
                                                  const DataToolTipsManager* dataToolTipsManager,
                                                  IdentificationStringBuilder& idText) const
{
    if (dataToolTipsManager->isShowChart()) {
        auto htmlTableBuilder = createHtmlTableBuilder(3);
        
        CaretMappableDataFile* dummyMapFile(NULL);
        std::set<int32_t> dummyMapIndices;
        this->generateChartTwoHistogramIdentificationText(*htmlTableBuilder,
                                                          idText,
                                                          selectionManager->getChartTwoHistogramIdentification(),
                                                          dummyMapFile,
                                                          dummyMapIndices,
                                                          true);

        this->generateChartTwoLineSeriesIdentificationText(*htmlTableBuilder,
                                                           idText,
                                                           selectionManager->getChartTwoLineSeriesIdentification(),
                                                           dummyMapFile,
                                                           dummyMapIndices,
                                                           true);

        this->generateChartTwoMatrixIdentificationText(*htmlTableBuilder,
                                                       idText,
                                                       selectionManager->getChartTwoMatrixIdentification(),
                                                       dummyMapFile,
                                                       dummyMapIndices,
                                                       true);
    }
}



/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
IdentificationFormattedTextGenerator::toString() const
{
    return "IdentificationFormattedTextGenerator";
}
        
/**
 * @return New instance of an HTML Table builder
 * @param numberOfColumns
 *  Number of columns in the table
 */
std::unique_ptr<HtmlTableBuilder>
IdentificationFormattedTextGenerator::createHtmlTableBuilder(const int32_t numberOfColumns) const
{
    std::unique_ptr<HtmlTableBuilder> htb(new HtmlTableBuilder(HtmlTableBuilder::HtmlVersion::V4_01,
                                                               numberOfColumns));
    return htb;
}
