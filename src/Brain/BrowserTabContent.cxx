
/*LICENSE_START*/
/* 
 *  Copyright 1995-2011 Washington University School of Medicine 
 * 
 *  http://brainmap.wustl.edu 
 * 
 *  This file is part of CARET. 
 * 
 *  CARET is free software; you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation; either version 2 of the License, or 
 *  (at your option) any later version. 
 * 
 *  CARET is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 * 
 *  You should have received a copy of the GNU General Public License 
 *  along with CARET; if not, write to the Free Software 
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
 * 
 */ 

#include <cmath>
#include <set>

#define __BROWSER_TAB_CONTENT_DECLARE__
#include "BrowserTabContent.h"
#undef __BROWSER_TAB_CONTENT_DECLARE__

#include "BorderFile.h"
#include "Brain.h"
#include "BrainOpenGLViewportContent.h"
#include "BrainStructure.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "GroupAndNameHierarchyGroup.h"
#include "GroupAndNameHierarchyModel.h"
#include "GroupAndNameHierarchyName.h"
#include "DisplayPropertiesBorders.h"
#include "DisplayPropertiesFoci.h"
#include "EventIdentificationHighlightLocation.h"
#include "EventModelGetAll.h"
#include "EventManager.h"
#include "FociFile.h"
#include "IdentificationManager.h"
#include "LabelFile.h"
#include "MathFunctions.h"
#include "Matrix4x4.h"
#include "ModelSurface.h"
#include "ModelSurfaceMontage.h"
#include "ModelSurfaceSelector.h"
#include "ModelTransform.h"
#include "ModelVolume.h"
#include "ModelWholeBrain.h"
#include "Overlay.h"
#include "OverlaySet.h"
#include "SceneAttributes.h"
#include "SceneClass.h"
#include "SceneClassAssistant.h"
#include "Surface.h"
#include "SurfaceSelectionModel.h"
#include "StructureEnum.h"
#include "VolumeFile.h"
#include "ViewingTransformations.h"
#include "VolumeSliceSettings.h"
#include "VolumeSurfaceOutlineModel.h"
#include "VolumeSurfaceOutlineSetModel.h"
#include "WholeBrainSurfaceSettings.h"

using namespace caret;

/**
 * Constructor.
 * @param tabNumber
 *    Number for this tab.
 */
BrowserTabContent::BrowserTabContent(const int32_t tabNumber)
: CaretObject()
{
    s_allBrowserTabContent.insert(this);
    
    m_tabNumber = tabNumber;
    m_surfaceModelSelector = new ModelSurfaceSelector();
    m_selectedModelType = ModelTypeEnum::MODEL_TYPE_INVALID;
    m_volumeModel = NULL;
    m_wholeBrainModel = NULL;
    m_surfaceMontageModel = NULL;
    m_guiName = "";
    m_userName = "";
    m_volumeSurfaceOutlineSetModel = new VolumeSurfaceOutlineSetModel();
    m_yokingGroup = YokingGroupEnum::YOKING_GROUP_OFF;

    m_viewingTransformation = new ViewingTransformations();
    m_volumeSliceViewingTransformation = new ViewingTransformations();
    m_wholeBrainSurfaceSettings = new WholeBrainSurfaceSettings();
    
    m_obliqueVolumeRotationMatrix = new Matrix4x4();
    
    leftView();

    m_volumeSliceSettings = new VolumeSliceSettings();
    
    m_clippingCoordinate[0] = 0.0;
    m_clippingCoordinate[1] = 0.0;
    m_clippingCoordinate[2] = 0.0;

    m_clippingThickness[0] = 20.0;
    m_clippingThickness[1] = 20.0;
    m_clippingThickness[2] = 20.0;
    
    m_clippingEnabled[0] = false;
    m_clippingEnabled[1] = false;
    m_clippingEnabled[2] = false;
    
    m_sceneClassAssistant = new SceneClassAssistant();
    m_sceneClassAssistant->add("m_tabNumber", 
                               &m_tabNumber);
    m_sceneClassAssistant->add("m_userName",
                               &m_userName);
    m_sceneClassAssistant->add<ModelTypeEnum, ModelTypeEnum::Enum>("m_selectedModelType",
                                                                   &m_selectedModelType);
    m_sceneClassAssistant->add("m_surfaceModelSelector",
                               "ModelSurfaceSelector",
                               m_surfaceModelSelector);
    m_sceneClassAssistant->add("m_volumeSurfaceOutlineSetModel",
                               "VolumeSurfaceOutlineSetModel",
                               m_volumeSurfaceOutlineSetModel);
    m_sceneClassAssistant->addArray("m_clippingCoordinate",
                                    m_clippingCoordinate,
                                    3,
                                    0.0);
    m_sceneClassAssistant->addArray("m_clippingThickness",
                                    m_clippingThickness,
                                    3,
                                    100.0);
    m_sceneClassAssistant->addArray("m_clippingEnabled",
                                    m_clippingEnabled,
                                    3,
                                    false);
    
    m_sceneClassAssistant->add("m_viewingTransformation",
                               "ViewingTransformations",
                               m_viewingTransformation);
    
    m_sceneClassAssistant->add("m_volumeSliceViewingTransformation",
                               "ViewingTransformations",
                               m_volumeSliceViewingTransformation);
    
    m_sceneClassAssistant->add("m_volumeSliceSettings",
                               "VolumeSliceSettings",
                               m_volumeSliceSettings);
    
    m_sceneClassAssistant->add("m_wholeBrainSurfaceSettings",
                               "WholeBrainSurfaceSettings",
                               m_wholeBrainSurfaceSettings);

    m_sceneClassAssistant->add<YokingGroupEnum, YokingGroupEnum::Enum>("m_yokingGroup",
                                                                   &m_yokingGroup);
    
    EventManager::get()->addEventListener(this,
                                          EventTypeEnum::EVENT_IDENTIFICATION_HIGHLIGHT_LOCATION);
}

/**
 * Destructor.
 */
BrowserTabContent::~BrowserTabContent()
{
    EventManager::get()->removeAllEventsFromListener(this);
 
    s_allBrowserTabContent.erase(this);
    
    delete m_viewingTransformation;
    delete m_volumeSliceViewingTransformation;
    delete m_obliqueVolumeRotationMatrix;
    
    delete m_surfaceModelSelector;
    m_surfaceModelSelector = NULL;
    
    delete m_volumeSurfaceOutlineSetModel;
    m_volumeSurfaceOutlineSetModel = NULL;
    
    delete m_volumeSliceSettings;
    
    delete m_wholeBrainSurfaceSettings;
    
    delete m_sceneClassAssistant;
    m_sceneClassAssistant = NULL;
}

/**
 * Clone the contents of the given browser tab.
 * @param tabToClone
 *    Tab whose contents is cloned.
 */
void 
BrowserTabContent::cloneBrowserTabContent(BrowserTabContent* tabToClone)
{
    CaretAssert(tabToClone);
    m_surfaceModelSelector->setSelectedStructure(tabToClone->m_surfaceModelSelector->getSelectedStructure());
    m_surfaceModelSelector->setSelectedSurfaceController(tabToClone->m_surfaceModelSelector->getSelectedSurfaceController());

    m_selectedModelType = tabToClone->m_selectedModelType;
    m_volumeModel = tabToClone->m_volumeModel;
    m_wholeBrainModel = tabToClone->m_wholeBrainModel;
    m_surfaceMontageModel = tabToClone->m_surfaceMontageModel;
    m_yokingGroup = tabToClone->m_yokingGroup;
    
    *m_viewingTransformation = *tabToClone->m_viewingTransformation;
    *m_volumeSliceViewingTransformation = *tabToClone->m_volumeSliceViewingTransformation;
    *m_volumeSliceSettings = *tabToClone->m_volumeSliceSettings;
    *m_wholeBrainSurfaceSettings = *tabToClone->m_wholeBrainSurfaceSettings;
    
    *m_obliqueVolumeRotationMatrix = *tabToClone->m_obliqueVolumeRotationMatrix;
    
    Model* model = getModelControllerForDisplay();
    
    const OverlaySet* overlaySetToClone = tabToClone->getOverlaySet();
    if (overlaySetToClone != NULL) {
        OverlaySet* overlaySet = getOverlaySet();
        if (overlaySet != NULL) {
            overlaySet->copyOverlaySet(overlaySetToClone);
        }
    }
    
    if (model != NULL) {
        Brain* brain = model->getBrain();
        
        brain->copyDisplayProperties(tabToClone->getTabNumber(),
                                     getTabNumber());
        
        const int32_t numberOfBrainStructures = brain->getNumberOfBrainStructures();
        for (int32_t i = 0; i < numberOfBrainStructures; i++) {
            BrainStructure* bs = brain->getBrainStructure(i);
            const int32_t numLabelFiles = bs->getNumberOfLabelFiles();
            for (int32_t j = 0; j < numLabelFiles; j++) {
                LabelFile* labelFile = bs->getLabelFile(j);
                labelFile->getGroupAndNameHierarchyModel()->copySelections(tabToClone->getTabNumber(),
                                                                           getTabNumber());
            }
        }

        const int32_t numBorderFiles = brain->getNumberOfBorderFiles();
        for (int32_t i = 0; i < numBorderFiles; i++) {
            BorderFile* bf = brain->getBorderFile(i);
            bf->getGroupAndNameHierarchyModel()->copySelections(tabToClone->getTabNumber(),
                                                                getTabNumber());
        }
        
        const int32_t numFociFiles = brain->getNumberOfFociFiles();
        for (int32_t i = 0; i < numFociFiles; i++) {
            FociFile* ff = brain->getFociFile(i);
            ff->getGroupAndNameHierarchyModel()->copySelections(tabToClone->getTabNumber(),
                                                                getTabNumber());
        }
    }
    
    m_volumeSurfaceOutlineSetModel->copyVolumeSurfaceOutlineSetModel(tabToClone->getVolumeSurfaceOutlineSet());
}

/**
 * Get the name of this browser tab.
 * Name priority is (1) name set by user, (2) name set by
 * user-interface, and (3) the default name.
 *
 * @return  Name of this tab.
 */
AString 
BrowserTabContent::getName() const
{
    AString s = "(" + AString::number(m_tabNumber + 1) + ") ";
    
    if (m_userName.isEmpty() == false) {
        s += m_userName;
    }
    else {
        const Model* displayedController =
            getModelControllerForDisplay();
        if (displayedController != NULL) {
            const AString name = displayedController->getNameForBrowserTab();
            s += name;
        }
    }
    
    return s;
}

/**
 * Set the user name of this tab.  The user name
 * overrides the default naming.
 *
 * @param userName
 *    User name for tab.
 */
void 
BrowserTabContent::setUserName(const AString& userName)
{
    m_userName = userName;
}

/**
 * @return The user name.
 */
AString 
BrowserTabContent::getUserName() const
{
    return m_userName;
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
BrowserTabContent::toString() const
{
    return "WindowTabContent";
}

/**
 * Get the selected model type.
 * 
 * @return The selected model type.
 */   
ModelTypeEnum::Enum 
BrowserTabContent::getSelectedModelType() const
{
    return m_selectedModelType;
}

/**
 * Set the selected model type.
 *
 * @param selectedModelType
 *    New selected model type. 
 */   
void 
BrowserTabContent::setSelectedModelType(ModelTypeEnum::Enum selectedModelType)
{
    m_selectedModelType = selectedModelType;
}

/**
 * Get the model controller for DISPLAY.
 * 
 * @return  Pointer to displayed controller or NULL
 *          if none are available.
 */   
Model* 
BrowserTabContent::getModelControllerForDisplay()
{
    Model* mdc = NULL;
    
    switch (m_selectedModelType) {
        case ModelTypeEnum::MODEL_TYPE_INVALID:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE:
            mdc = m_surfaceModelSelector->getSelectedSurfaceController();
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE:
            mdc = m_surfaceMontageModel;
            break;
        case ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES:
            mdc = m_volumeModel;
            break;
        case ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN:
            mdc = m_wholeBrainModel;
            break;
        case ModelTypeEnum::MODEL_TYPE_CHART:
            break;
    }
    
    return mdc;
}

/**
 * Get the model controller for DISPLAY.
 * 
 * @return  Pointer to displayed controller or NULL
 *          if none are available.
 */   
const Model* 
BrowserTabContent::getModelControllerForDisplay() const
{
    Model* mdc = NULL;
    
    switch (m_selectedModelType) {
        case ModelTypeEnum::MODEL_TYPE_INVALID:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE:
            mdc = m_surfaceModelSelector->getSelectedSurfaceController();
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE:
            mdc = m_surfaceMontageModel;
            break;
        case ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES:
            mdc = m_volumeModel;
            break;
        case ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN:
            mdc = m_wholeBrainModel;
            break;
        case ModelTypeEnum::MODEL_TYPE_CHART:
            break;
    }
    
    return mdc;
}

/**
 * Get the displayed surface model.
 * 
 * @return  Pointer to displayed surface model or 
 *          NULL if the displayed model is NOT a 
 *          surface.
 */   
ModelSurface* 
BrowserTabContent::getDisplayedSurfaceModel()
{
    ModelSurface* mdcs =
        dynamic_cast<ModelSurface*>(getModelControllerForDisplay());
    return mdcs;
}

/**
 * Get the displayed surface model.
 * 
 * @return  Pointer to displayed surface model or 
 *          NULL if the displayed model is NOT a 
 *          surface.
 */   
const ModelSurface* 
BrowserTabContent::getDisplayedSurfaceModel() const
{
    const ModelSurface* mdcs =
    dynamic_cast<const ModelSurface*>(getModelControllerForDisplay());
    return mdcs;
}

/**
 * Get the displayed volume model.
 * 
 * @return  Pointer to displayed volume model or 
 *          NULL if the displayed model is NOT a 
 *          volume.
 */   
ModelVolume* 
BrowserTabContent::getDisplayedVolumeModel()
{
    ModelVolume* mdcv =
        dynamic_cast<ModelVolume*>(getModelControllerForDisplay());
    return mdcv;
}

/**
 * @return Is the displayed model a volume slice model?
 */
bool
BrowserTabContent::isVolumeSlicesDisplayed() const
{
    const ModelVolume* mdcv = dynamic_cast<const ModelVolume*>(getModelControllerForDisplay());
    
    const bool volumeFlag = (mdcv != NULL);
    return volumeFlag;
}

/**
 * @return Is the displayed model the whole brain model (ALL)?
 */
bool
BrowserTabContent::isWholeBrainDisplayed() const
{
    const ModelWholeBrain* mwb = dynamic_cast<const ModelWholeBrain*>(getModelControllerForDisplay());
    const bool wholeBrainFlag = (mwb != NULL);
    return wholeBrainFlag;
}

/**
 * Get the displayed whole brain model.
 * 
 * @return  Pointer to displayed whole brain model or 
 *          NULL if the displayed model is NOT a 
 *          whole brain.
 */   
ModelWholeBrain* 
BrowserTabContent::getDisplayedWholeBrainModel()
{
    ModelWholeBrain* mdcwb =
        dynamic_cast<ModelWholeBrain*>(getModelControllerForDisplay());
    return mdcwb;

}

/**
 * @return Pointer to displayed surface montage model
 * or NULL if the displayed model is not a surface
 * montage model.
 */
ModelSurfaceMontage* 
BrowserTabContent::getDisplayedSurfaceMontageModel()
{
    ModelSurfaceMontage* mdcsm =
    dynamic_cast<ModelSurfaceMontage*>(getModelControllerForDisplay());
    return mdcsm;
}

/**
 * Get all of the available surface models.
 * 
 * @return Vector containing all surface models.
 */   
const std::vector<ModelSurface*> 
BrowserTabContent::getAllSurfaceModels() const
{
    return m_allSurfaceModels;
}

/**
 * @return The surface model selector used to 
 * select a surface and structure.
 */
ModelSurfaceSelector* 
BrowserTabContent::getSurfaceModelSelector()
{
    return m_surfaceModelSelector;
}

/**
 * Get the overlay assignments for this tab.
 * 
 * @return  Overlay assignments for this tab or NULL if no valid model.
 */
OverlaySet* 
BrowserTabContent::getOverlaySet()
{
    Model* modelDisplayController = getModelControllerForDisplay();
    if (modelDisplayController != NULL) {
       return modelDisplayController->getOverlaySet(m_tabNumber);
    }
    return NULL;
}

/**
 * Get the tab number for this content.
 * 
 * @return  Tab number.
 */
int32_t 
BrowserTabContent::getTabNumber() const
{
    return m_tabNumber;
}

/**
 * Update the selected models.
 */
void 
BrowserTabContent::update(const std::vector<Model*> modelDisplayControllers)
{
    m_surfaceModelSelector->updateSelector(modelDisplayControllers);
    
    const int32_t numModels = static_cast<int32_t>(modelDisplayControllers.size());
    
    ModelVolume* previousVolumeModel = m_volumeModel;
    
    m_allSurfaceModels.clear();
    m_surfaceModelSelector->getSelectableSurfaceControllers(m_allSurfaceModels);
    m_volumeModel = NULL;
    m_wholeBrainModel = NULL;
    m_surfaceMontageModel = NULL;
    
    for (int i = 0; i < numModels; i++) {
        Model* mdc = modelDisplayControllers[i];
        
        ModelSurface* mdcs = dynamic_cast<ModelSurface*>(mdc);
        ModelVolume* mdcv = dynamic_cast<ModelVolume*>(mdc);
        ModelWholeBrain* mdcwb = dynamic_cast<ModelWholeBrain*>(mdc);
        ModelSurfaceMontage* mdcsm = dynamic_cast<ModelSurfaceMontage*>(mdc);
        
        if (mdcs != NULL) {
            /* nothing to do since the surface model selector handles surfaces */
        }
        else if (mdcv != NULL) {
            CaretAssertMessage((m_volumeModel == NULL), "There is more than one volume model.");
            m_volumeModel = mdcv;
        }
        else if (mdcwb != NULL) {
            CaretAssertMessage((m_wholeBrainModel == NULL), "There is more than one whole brain model.");
            m_wholeBrainModel = mdcwb;
        }
        else if (mdcsm != NULL) {
            CaretAssertMessage((m_surfaceMontageModel == NULL), "There is more than one surface montage model.");
            m_surfaceMontageModel = mdcsm;
        }
        else {
            CaretAssertMessage(0, (AString("Unknown type of brain model.") + mdc->getNameForGUI(true)));
        }
    }
    
    switch (m_selectedModelType) {
        case ModelTypeEnum::MODEL_TYPE_INVALID:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE:
            if (m_surfaceModelSelector->getSelectedSurfaceController() == NULL) {
                m_selectedModelType = ModelTypeEnum::MODEL_TYPE_INVALID;
            }
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE:
            if (m_surfaceMontageModel == NULL) {
                m_selectedModelType = ModelTypeEnum::MODEL_TYPE_INVALID;
            }
            break;
        case ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES:
            if (m_volumeModel == NULL) {
                m_selectedModelType = ModelTypeEnum::MODEL_TYPE_INVALID;
            }
            break;
        case ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN:
            if (m_wholeBrainModel == NULL) {
                m_selectedModelType = ModelTypeEnum::MODEL_TYPE_INVALID;
            }
            break;
        case ModelTypeEnum::MODEL_TYPE_CHART:
            break;
    }
    
    if (m_selectedModelType == ModelTypeEnum::MODEL_TYPE_INVALID) {
        if (m_surfaceModelSelector->getSelectedSurfaceController() != NULL) {
            m_selectedModelType = ModelTypeEnum::MODEL_TYPE_SURFACE;
        }
        else if (m_volumeModel != NULL) {
            m_selectedModelType = ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES;
        }
        else if (m_wholeBrainModel != NULL) {
            m_selectedModelType = ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN;
        }
        else if (m_surfaceMontageModel != NULL) {
            m_selectedModelType = ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE;
        }
    }
    
    if (m_volumeModel != NULL) {
        if (m_volumeModel != previousVolumeModel) {
        
            VolumeMappableInterface* underlayVolume = m_volumeModel->getOverlaySet(m_tabNumber)->getUnderlayVolume();
            if (underlayVolume != NULL) {
                /*
                 * Set montage slice spacing based upon slices
                 * in the Z-axis.
                 */
                std::vector<int64_t> dimensions;
                underlayVolume->getDimensions(dimensions);
                
                if (dimensions.size() >= 3) {
                    const int32_t dimZ = dimensions[2];
                    if (dimZ > 0) {
                        const int32_t maxZ = static_cast<int32_t>(dimZ * 0.90);
                        const int32_t minZ = static_cast<int32_t>(dimZ * 0.10);
                        const int32_t sliceRange = (maxZ - minZ);
                        int32_t sliceSpacing = 1;
                        if (sliceRange > 0) {
                            const int32_t numSlicesViewed = (m_volumeSliceSettings->getMontageNumberOfRows()
                                                             * m_volumeSliceSettings->getMontageNumberOfColumns());
                            sliceSpacing = (sliceRange / numSlicesViewed);
                        }
                        m_volumeSliceSettings->setMontageSliceSpacing(sliceSpacing);
                    }
                }
            }
        }
    }
}

/**
 * Is the surface model selection valid?
 *
 * @return bool indicating validity.
 */
bool 
BrowserTabContent::isSurfaceModelValid() const
{
    bool valid = (m_allSurfaceModels.empty() == false);
    return valid;
}

/**
 * Is the volume model selection valid?
 *
 * @return bool indicating validity.
 */
bool 
BrowserTabContent::isVolumeSliceModelValid() const
{
    bool valid = (m_volumeModel != NULL);
    return valid;
}

/**
 * Is the whole brain model selection valid?
 *
 * @return bool indicating validity.
 */
bool 
BrowserTabContent::isWholeBrainModelValid() const
{
    bool valid = (m_wholeBrainModel != NULL);
    return valid;
}

/**
 * Is the surface montage model selection valid?
 *
 * @return bool indicating validity.
 */
bool 
BrowserTabContent::isSurfaceMontageModelValid() const
{
    bool valid = (m_surfaceMontageModel != NULL);
    return valid;
}

/**
 * Receive an event.
 * 
 * @param event
 *     The event that the receive can respond to.
 */
void
BrowserTabContent::receiveEvent(Event* event)
{
    if (event->getEventType() == EventTypeEnum::EVENT_IDENTIFICATION_HIGHLIGHT_LOCATION) {
        EventIdentificationHighlightLocation* idLocationEvent =
        dynamic_cast<EventIdentificationHighlightLocation*>(event);
        CaretAssert(idLocationEvent);

        Model* model = getModelControllerForDisplay();
        if (model == NULL) {
            return;
        }

        Brain* brain = model->getBrain();
        if (brain->getIdentificationManager()->isVolumeIdentificationEnabled()) {
            const float* highlighXYZ = idLocationEvent->getXYZ();
            for (int32_t windowTabNumber = 0;
                 windowTabNumber < BrainConstants::MAXIMUM_NUMBER_OF_BROWSER_TABS;
                 windowTabNumber++) {
                
                float volumeSliceXYZ[3] = {
                    highlighXYZ[0],
                    highlighXYZ[1],
                    highlighXYZ[2]
                };
                
                /*
                 * If volume montage viewing, do not change the
                 * slice in the plane that is being viewed.
                 */
                if (getDisplayedVolumeModel() != NULL) {
                    switch (m_volumeSliceSettings->getSliceViewMode()) {
                        case VolumeSliceViewModeEnum::MONTAGE:
                            switch (getSliceViewPlane()) {
                                case VolumeSliceViewPlaneEnum::ALL:
                                    break;
                                case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                                    volumeSliceXYZ[0] = getSliceCoordinateParasagittal();
                                    break;
                                case VolumeSliceViewPlaneEnum::CORONAL:
                                    volumeSliceXYZ[1] = getSliceCoordinateCoronal();
                                    break;
                                case VolumeSliceViewPlaneEnum::AXIAL:
                                    volumeSliceXYZ[2] = getSliceCoordinateAxial();
                                    break;
                            }
                            break;
                        case VolumeSliceViewModeEnum::OBLIQUE:
                            break;
                        case VolumeSliceViewModeEnum::ORTHOGONAL:
                            break;
                    }
                }
                
                m_volumeSliceSettings->selectSlicesAtCoordinate(volumeSliceXYZ);
            }
        }
        
        idLocationEvent->setEventProcessed();
    }
}

/**
 * Get the map files for which a palette should be displayed in the 
 * graphcis window.  Note that the order of map files and indices starts
 * with the bottom most layer and ends with the top most overlay.
 *
 * @param mapFiles
 *    Outut Map files that should have a palette displayed in the graphics window.
 * @param mapIndices
 *    Output Indices of the maps in the mapFiles.
 */
void 
BrowserTabContent::getDisplayedPaletteMapFiles(std::vector<CaretMappableDataFile*>& mapFiles,
                                               std::vector<int32_t>& mapIndices)
{
    mapFiles.clear();
    mapIndices.clear();
    
    OverlaySet* overlaySet = getOverlaySet();
    const int32_t numOverlays = overlaySet->getNumberOfDisplayedOverlays();
    for (int32_t i = (numOverlays - 1); i >= 0; i--) {
        Overlay* overlay = overlaySet->getOverlay(i);
        if (overlay->isEnabled()) {
            if (overlay->isPaletteDisplayEnabled()) {
                CaretMappableDataFile* mapFile;
                int32_t mapFileIndex;
                overlay->getSelectionData(mapFile, mapFileIndex);
                if (mapFile != NULL) {
                    if (mapFile->isMappedWithPalette()) {
                        if ((mapFileIndex >= 0) 
                            && (mapFileIndex < mapFile->getNumberOfMaps())) {
                            mapFiles.push_back(mapFile);
                            mapIndices.push_back(mapFileIndex);
                        }
                    }
                }
            }
        }
    }
}

/**
 * @return The volume surface outline model for this tab.
 */
VolumeSurfaceOutlineSetModel* 
BrowserTabContent::getVolumeSurfaceOutlineSet()
{
    return m_volumeSurfaceOutlineSetModel;
}

/**
 * @return The volume surface outline model for this tab.
 */
const VolumeSurfaceOutlineSetModel* 
BrowserTabContent::getVolumeSurfaceOutlineSet() const
{
    return m_volumeSurfaceOutlineSetModel;
}

/**
 * Get the data files displayed in this tab.
 * @param displayedDataFilesOut
 *    Displayed data file are loaded into this parameter.
 */
void 
BrowserTabContent::getFilesDisplayedInTab(std::vector<CaretDataFile*>& displayedDataFilesOut)
{
    displayedDataFilesOut.clear();
 
    Model* model = getModelControllerForDisplay();
    if (model == NULL) {
        return;
    }
    
    std::set<CaretDataFile*> displayedDataFiles;
    
    const int32_t tabIndex = getTabNumber();
    switch (getSelectedModelType()) {
        case ModelTypeEnum::MODEL_TYPE_INVALID:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE:
        {
            ModelSurface* ms = getDisplayedSurfaceModel();
            displayedDataFiles.insert(ms->getSurface());
        }
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE:
        {
            ModelSurfaceMontage* msm = getDisplayedSurfaceMontageModel();
            if (msm->isFirstSurfaceEnabled(tabIndex)) {
                if (msm->isLeftEnabled(tabIndex)) {
                    displayedDataFiles.insert(msm->getLeftSurfaceSelectionModel(tabIndex)->getSurface());
                }
                if (msm->isRightEnabled(tabIndex)) {
                    displayedDataFiles.insert(msm->getRightSurfaceSelectionModel(tabIndex)->getSurface());
                }
            }
            if (msm->isSecondSurfaceEnabled(tabIndex)) {
                if (msm->isLeftEnabled(tabIndex)) {
                    displayedDataFiles.insert(msm->getLeftSecondSurfaceSelectionModel(tabIndex)->getSurface());
                }
                if (msm->isRightEnabled(tabIndex)) {
                    displayedDataFiles.insert(msm->getRightSecondSurfaceSelectionModel(tabIndex)->getSurface());
                }
            }
        }
            break;
        case ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES:
        {
            const int32_t numOutlines = m_volumeSurfaceOutlineSetModel->getNumberOfDislayedVolumeSurfaceOutlines();
            for (int32_t i = 0; i < numOutlines; i++) {
                VolumeSurfaceOutlineModel* model = m_volumeSurfaceOutlineSetModel->getVolumeSurfaceOutlineModel(i);
                if (model->isDisplayed()) {
                    displayedDataFiles.insert(model->getSurface());
                }
            }
        }
            break;
        case ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN:
        {
            ModelWholeBrain* wbm = getDisplayedWholeBrainModel();
            if (isWholeBrainLeftEnabled()) {
                displayedDataFiles.insert(wbm->getSelectedSurface(StructureEnum::CORTEX_LEFT, tabIndex));
            }
            if (isWholeBrainRightEnabled()) {
                displayedDataFiles.insert(wbm->getSelectedSurface(StructureEnum::CORTEX_RIGHT, tabIndex));
            }
            if (isWholeBrainCerebellumEnabled()) {
                displayedDataFiles.insert(wbm->getSelectedSurface(StructureEnum::CEREBELLUM, tabIndex));
            }
        }
            break;
        case ModelTypeEnum::MODEL_TYPE_CHART:
            break;
    }

    /*
     * Check overlay files
     */
    OverlaySet* overlaySet = model->getOverlaySet(tabIndex);
    const int32_t numOverlays = overlaySet->getNumberOfDisplayedOverlays();
    for (int32_t i = 0; i < numOverlays; i++) {
        Overlay* overlay = overlaySet->getOverlay(i);
        CaretMappableDataFile* overlayDataFile = NULL;
        int32_t mapIndex;
        overlay->getSelectionData(overlayDataFile, 
                                  mapIndex);
        displayedDataFiles.insert(overlayDataFile);
    }
    
    
    /*
     * Check border files
     */
    Brain* brain = model->getBrain();
    const DisplayPropertiesBorders* borderDisplayProperties = brain->getDisplayPropertiesBorders();
    const DisplayGroupEnum::Enum borderDisplayGroup = borderDisplayProperties->getDisplayGroupForTab(tabIndex);
    if (borderDisplayProperties->isDisplayed(borderDisplayGroup,
                                             tabIndex)) {
        const int32_t numBorderFiles = brain->getNumberOfBorderFiles();
        for (int32_t i = 0; i < numBorderFiles; i++) {
            BorderFile* borderFile = brain->getBorderFile(i);
            const GroupAndNameHierarchyModel* classAndNameSelection = borderFile->getGroupAndNameHierarchyModel();
            if (classAndNameSelection->isSelected(borderDisplayGroup,
                                                  tabIndex)) {
                displayedDataFilesOut.push_back(borderFile);
            }
        }
    }

    /*
     * Check foci files
     */
    const DisplayPropertiesFoci* fociDisplayProperties = brain->getDisplayPropertiesFoci();
    const DisplayGroupEnum::Enum fociDisplayGroup = fociDisplayProperties->getDisplayGroupForTab(tabIndex);
    
    if (fociDisplayProperties->isDisplayed(fociDisplayGroup,
                                           tabIndex)) {
        const int32_t numFociFiles = brain->getNumberOfFociFiles();
        for (int32_t i = 0; i < numFociFiles; i++) {
            FociFile* fociFile = brain->getFociFile(i);
            
            const GroupAndNameHierarchyModel* classAndNameSelection = fociFile->getGroupAndNameHierarchyModel();
            if (classAndNameSelection->isSelected(fociDisplayGroup,
                                                  tabIndex)) {
                displayedDataFilesOut.push_back(fociFile);
            }
        }
    }
    
    /*
     * Might be NULLs so filter them out and return the results
     */
    for (std::set<CaretDataFile*>::iterator iter = displayedDataFiles.begin();
         iter != displayedDataFiles.end();
         iter++) {
        CaretDataFile* cdf = *iter;
        if (cdf != NULL) {
            displayedDataFilesOut.push_back(cdf);
        }
    }
}

/**
 * @return The viewing translation.
 */
const float*
BrowserTabContent::getTranslation() const
{
    if (isVolumeSlicesDisplayed()) {
        return m_volumeSliceViewingTransformation->getTranslation();
    }
    else {
        return m_viewingTransformation->getTranslation();
    }
}

/**
 * Get the viewing translation.
 *
 * @param translationOut
 *    Translation values output.
 */
void
BrowserTabContent::getTranslation(float translationOut[3]) const
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->getTranslation(translationOut);
    }
    else {
        m_viewingTransformation->getTranslation(translationOut);
    }
}

/**
 * Set the viewing translation.
 *
 * @param translation
 *    New translation values.
 */
void
BrowserTabContent::setTranslation( const float translation[3])
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->setTranslation(translation);
    }
    else {
        m_viewingTransformation->setTranslation(translation);
    }
    updateYokedBrowserTabs();
}

/**
 * Set the viewing translation.
 *
 * @param translationX
 *    New translation X-value.
 * @param translationY
 *    New translation Y-value.
 * @param translationZ
 *    New translation Z-value.
 */
void
BrowserTabContent::setTranslation(const float translationX,
                                  const float translationY,
                                  const float translationZ)
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->setTranslation(translationX,
                                                           translationY,
                                                           translationZ);
    }
    else {
        m_viewingTransformation->setTranslation(translationX,
                                                translationY,
                                                translationZ);
    }
    updateYokedBrowserTabs();
}

/**
 * @return The viewing scaling.
 */
float
BrowserTabContent::getScaling() const
{
    if (isVolumeSlicesDisplayed()) {
        return m_volumeSliceViewingTransformation->getScaling();
    }
    else {
        return m_viewingTransformation->getScaling();
    }
}

/**
 * Set the viewing scaling.
 * @param scaling
 *    New value for scaling.
 */
void
BrowserTabContent::setScaling(const float scaling)
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->setScaling(scaling);
    }
    else {
        m_viewingTransformation->setScaling(scaling);
    }
    updateYokedBrowserTabs();
}

/**
 * @return The rotation matrix.
 */
Matrix4x4
BrowserTabContent::getRotationMatrix() const
{
    if (isVolumeSlicesDisplayed()) {
        return m_volumeSliceViewingTransformation->getRotationMatrix();
    }
    else {
        return m_viewingTransformation->getRotationMatrix();
    }
}

/**
 * Set the rotation matrix.
 *
 * @param rotationMatrix
 *    The new rotation matrix.
 */
void
BrowserTabContent::setRotationMatrix(const Matrix4x4& rotationMatrix)
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->setRotationMatrix(rotationMatrix);
    }
    else {
        m_viewingTransformation->setRotationMatrix(rotationMatrix);
    }
    updateYokedBrowserTabs();
}

/**
 * @return The oblique volume rotation matrix.
 */
Matrix4x4
BrowserTabContent::getObliqueVolumeRotationMatrix() const
{
    return *m_obliqueVolumeRotationMatrix;
}

/**
 * Set the oblique rotation matrix.
 *
 * @param obliqueRotationMatrix
 *    The new oblique rotation matrix.
 */
void
BrowserTabContent::setObliqueVolumeRotationMatrix(const Matrix4x4& obliqueRotationMatrix)
{
    *m_obliqueVolumeRotationMatrix = obliqueRotationMatrix;
    
    updateYokedBrowserTabs();
}

/**
 * Reset the view to the default view.
 */
void
BrowserTabContent::resetView()
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->resetVolumeView();
        m_obliqueVolumeRotationMatrix->identity();
    }
    else {
        m_viewingTransformation->resetView();
    }
    updateYokedBrowserTabs();
}

/**
 * Set to a right side view.
 */
void
BrowserTabContent::rightView()
{
    if (isVolumeSlicesDisplayed()) {
        /* Nothing */
    }
    else {
        m_viewingTransformation->rightView();
    }
    updateYokedBrowserTabs();
}

/**
 * set to a left side view.
 */
void
BrowserTabContent::leftView()
{
    if (isVolumeSlicesDisplayed()) {
        /* Nothing */
    }
    else {
        m_viewingTransformation->leftView();
    }
    updateYokedBrowserTabs();
}

/**
 * set to a anterior view.
 */
void
BrowserTabContent::anteriorView()
{
    if (isVolumeSlicesDisplayed()) {
        /* Nothing */
    }
    else {
        m_viewingTransformation->anteriorView();
    }
    updateYokedBrowserTabs();
}

/**
 * set to a posterior view.
 */
void
BrowserTabContent::posteriorView()
{
    if (isVolumeSlicesDisplayed()) {
        /* Nothing */
    }
    else {
        m_viewingTransformation->posteriorView();
    }
    updateYokedBrowserTabs();
}

/**
 * set to a dorsal view.
 */
void
BrowserTabContent::dorsalView()
{
    if (isVolumeSlicesDisplayed()) {
        /* Nothing */
    }
    else {
        m_viewingTransformation->dorsalView();
    }
    updateYokedBrowserTabs();
}

/**
 * set to a ventral view.
 */
void
BrowserTabContent::ventralView()
{
    if (isVolumeSlicesDisplayed()) {
        /* Nothing */
    }
    else {
        m_viewingTransformation->ventralView();
    }
    updateYokedBrowserTabs();
}

/*
 * @return The slice view plane for the given viewport coordinate.
 * If ALL is returned, is indicates that the given viewport coordinate
 * is in the bottom left region in which volume slices are not displayed.
 *
 * @param viewport
 *   The viewport.
 * @param mousePressX
 *   X Location of the mouse press.
 * @param mousePressY
 *   Y Location of the mouse press.
 */
VolumeSliceViewPlaneEnum::Enum
BrowserTabContent::getSliceViewPlaneForVolumeAllSliceView(const int32_t viewport[4],
                                                          const int32_t mousePressX,
                                                          const int32_t mousePressY,
                                                          int32_t sliceViewportOut[4]) const
{
    VolumeSliceViewPlaneEnum::Enum view = VolumeSliceViewPlaneEnum::ALL;
    
    const int32_t halfWidth  = viewport[2] / 2;
    const int32_t halfHeight = viewport[3] / 2;
    const int32_t viewportMousePressedX = mousePressX - viewport[0];
    const int32_t viewportMousePressedY = mousePressY - viewport[1];
    bool isRight  = false;
    bool isTop = false;
    if (viewportMousePressedX > halfWidth) {
        isRight = true;
    }
    if (viewportMousePressedY > halfHeight) {
        isTop = true;
    }
    
    /*
     * Top Right is Coronal
     * Top Left is Parasagittal
     * Bottom Right is Axial
     * Bottom Left is Empty or Surfaces
     */
    if (isTop) {
        if (isRight) {
            view = VolumeSliceViewPlaneEnum::CORONAL;
        }
        else {
            view = VolumeSliceViewPlaneEnum::PARASAGITTAL;
        }
    }
    else {
        if (isRight) {
            view = VolumeSliceViewPlaneEnum::AXIAL;
        }
        else {
            
        }
    }

    sliceViewportOut[0] = viewport[0];
    sliceViewportOut[1] = viewport[1];
    sliceViewportOut[2] = halfWidth;
    sliceViewportOut[3] = halfHeight;
    
    switch (view) {
        case VolumeSliceViewPlaneEnum::ALL:
            sliceViewportOut[2] = viewport[2];
            sliceViewportOut[3] = viewport[3];
            break;
        case VolumeSliceViewPlaneEnum::AXIAL:
            sliceViewportOut[0] = halfWidth;
            break;
        case VolumeSliceViewPlaneEnum::CORONAL:
            sliceViewportOut[0] = halfWidth;
            sliceViewportOut[1] = halfHeight;
            break;
        case VolumeSliceViewPlaneEnum::PARASAGITTAL:
            sliceViewportOut[1] = halfHeight;
            break;
    }
    return view;
}

/**
 * Apply mouse rotation to the displayed model.
 *
 * @param mousePressX
 *    X coordinate of where mouse was pressed.
 * @param mousePressY
 *    X coordinate of where mouse was pressed.
 * @param mouseX
 *    X coordinate of mouse.
 * @param mouseY
 *    Y coordinate of mouse.
 * @param mouseDeltaX
 *    Change in mouse X coordinate.
 * @param mouseDeltaY
 *    Change in mouse Y coordinate.
 */
void
BrowserTabContent::applyMouseRotation(BrainOpenGLViewportContent* viewportContent,
                                      const int32_t mousePressX,
                                      const int32_t mousePressY,
                                      const int32_t mouseX,
                                      const int32_t mouseY,
                                      const int32_t mouseDeltaX,
                                      const int32_t mouseDeltaY)
{
    if (isVolumeSlicesDisplayed()) {
        switch (getSliceViewMode()) {
            case VolumeSliceViewModeEnum::MONTAGE:
                break;
            case VolumeSliceViewModeEnum::OBLIQUE:
            {
                int viewport[4];
                viewportContent->getModelViewport(viewport);
                VolumeSliceViewPlaneEnum::Enum slicePlane = this->getSliceViewPlane();
                int sliceViewport[4] = {
                    viewport[0],
                    viewport[1],
                    viewport[2],
                    viewport[3]
                };
                if (slicePlane == VolumeSliceViewPlaneEnum::ALL) {
                    slicePlane = getSliceViewPlaneForVolumeAllSliceView(viewport,
                                                                        mousePressX,
                                                                        mousePressY,
                                                                        sliceViewport);
                }
                
                Matrix4x4 rotationMatrix = getObliqueVolumeRotationMatrix(); //  m_volumeSliceViewingTransformation->getRotationMatrix();
                
                if (slicePlane == VolumeSliceViewPlaneEnum::ALL) {
                    rotationMatrix.rotateX(-mouseDeltaY);
                    rotationMatrix.rotateY(mouseDeltaX);
                }
                else {
                    if ((mouseDeltaX != 0)
                        && (mouseDeltaY != 0)) {
                        
                        const int previousMouseX = mouseX - mouseDeltaX;
                        const int previousMouseY = mouseY - mouseDeltaY;
                        
                        /*
                         * Need to account for the quadrants!!!!
                         */
                        const float viewportCenter[3] = {
                            sliceViewport[0] + sliceViewport[2] / 2,
                            sliceViewport[1] + sliceViewport[3] / 2,
                            0.0
                        };
                        
                        const float oldPos[3] = {
                            previousMouseX,
                            previousMouseY,
                            0.0
                        };
                        
                        const float newPos[3] = {
                            mouseX,
                            mouseY,
                            0.0
                        };
                        
//                        const AString posString = (" center: "
//                                             + AString::fromNumbers(viewportCenter, 3, ",")
//                                             + " old: "
//                                             + AString::fromNumbers(oldPos, 3, ",")
//                                             + " mouse: "
//                                             + AString::fromNumbers(newPos, 3, ","));
                        /*
                         * Compute normal vector from viewport center to
                         * old mouse position to new mouse position.
                         * If normal-Z is positive, mouse has been moved
                         * in a counter clockwise motion relative to center.
                         * If normal-Z is negative, mouse has moved clockwise.
                         */
                        float normalDirection[3];
                        MathFunctions::normalVectorDirection(viewportCenter,
                                                             oldPos,
                                                             newPos,
                                                             normalDirection);
                        bool isClockwise = false;
                        bool isCounterClockwise = false;
                        if (normalDirection[2] > 0.0) {
                            isCounterClockwise = true;
//                            std::cout << "CCW" << qPrintable(posString) << std::endl;
                        }
                        else if (normalDirection[2] < 0.0) {
                            isClockwise = true;
//                            std::cout << "CW" << qPrintable(posString) << std::endl;
                        }
                        
                        if (isClockwise
                            || isCounterClockwise) {
                            float mouseDelta = std::sqrt((mouseDeltaX * mouseDeltaX)
                                                         + (mouseDeltaY * mouseDeltaY));
                            
                            switch (slicePlane) {
                                case VolumeSliceViewPlaneEnum::ALL:
                                {
                                    CaretAssert(0);
                                }
                                    break;
                                case VolumeSliceViewPlaneEnum::AXIAL:
                                {
                                    Matrix4x4 rotation;
                                    if (isClockwise) {
                                        rotation.rotateZ(mouseDelta);
                                    }
                                    else if (isCounterClockwise) {
                                        rotation.rotateZ(-mouseDelta);
                                    }
                                    rotationMatrix.premultiply(rotation);
                                }
                                    break;
                                case VolumeSliceViewPlaneEnum::CORONAL:
                                {
                                    Matrix4x4 rotation;
                                    if (isClockwise) {
                                        rotation.rotateY(-mouseDelta);
                                    }
                                    else if (isCounterClockwise) {
                                        rotation.rotateY(mouseDelta);
                                    }
                                    rotationMatrix.premultiply(rotation);
                                }
                                    break;
                                case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                                {
                                    Matrix4x4 rotation;
                                    if (isClockwise) {
                                        rotation.rotateX(-mouseDelta);
                                    }
                                    else if (isCounterClockwise) {
                                        rotation.rotateX(mouseDelta);
                                    }
                                    rotationMatrix.premultiply(rotation);
                                }
                                    break;
                            }
                        }
                    }
                }
                
                setObliqueVolumeRotationMatrix(rotationMatrix);
                
//                if (mouseDeltaY != 0) {
//                    Matrix4x4 rotationMatrix = getObliqueVolumeRotationMatrix(); //  m_volumeSliceViewingTransformation->getRotationMatrix();
//                    switch (slicePlane) {
//                        case VolumeSliceViewPlaneEnum::ALL:
//                        {
//                            rotationMatrix.rotateX(-mouseDeltaY);
//                            rotationMatrix.rotateY(mouseDeltaX);
//                        }
//                            break;
//                        case VolumeSliceViewPlaneEnum::AXIAL:
//                        {
//                            Matrix4x4 rotation;
//                            rotation.rotateZ(mouseDeltaY);
//                            rotationMatrix.premultiply(rotation);
//                        }
//                            break;
//                        case VolumeSliceViewPlaneEnum::CORONAL:
//                        {
//                            Matrix4x4 rotation;
//                            rotation.rotateY(mouseDeltaY);
//                            rotationMatrix.premultiply(rotation);
//                        }
//                            break;
//                        case VolumeSliceViewPlaneEnum::PARASAGITTAL:
//                        {
//                            Matrix4x4 rotation;
//                            rotation.rotateX(mouseDeltaY);
//                            rotationMatrix.premultiply(rotation);
//                        }
//                            break;
//                    }
//                    setObliqueVolumeRotationMatrix(rotationMatrix);
////                    m_volumeSliceViewingTransformation->setRotationMatrix(rotationMatrix);
//                    
////                    Matrix4x4 rotation;
////                    switch (slicePlane) {
////                        case VolumeSliceViewPlaneEnum::ALL:
////                            rotation.rotateX(-mouseDY);
////                            rotation.rotateY(mouseDX);
////                            break;
////                        case VolumeSliceViewPlaneEnum::AXIAL:
////                            rotation.rotateZ(mouseDY);
////                            break;
////                        case VolumeSliceViewPlaneEnum::CORONAL:
////                            rotation.rotateY(mouseDY);
////                            break;
////                        case VolumeSliceViewPlaneEnum::PARASAGITTAL:
////                            rotation.rotateX(mouseDY);
////                            break;
////                    }
////
////                    Matrix4x4 rotationMatrix = m_volumeSliceViewingTransformation->getRotationMatrix();
////                    rotationMatrix.premultiply(rotation);
////                    m_volumeSliceViewingTransformation->setRotationMatrix(rotationMatrix);
//                }
                
            }
                break;
            case VolumeSliceViewModeEnum::ORTHOGONAL:
                /* Orthogonal olume slices are not rotated */
                break;
        }
    }
    else {
        if (getProjectionViewType() == ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_LATERAL) {
            Matrix4x4 rotationMatrix = m_viewingTransformation->getRotationMatrix();
            rotationMatrix.rotateX(-mouseDeltaY);
            rotationMatrix.rotateY(-mouseDeltaX);
            m_viewingTransformation->setRotationMatrix(rotationMatrix);
        }
        else {
            float dx = mouseDeltaX;
            float dy = mouseDeltaY;
            
            ModelSurfaceMontage* montageModel = getDisplayedSurfaceMontageModel();
            if (montageModel != NULL) {
                std::vector<SurfaceMontageViewport> montageViewports;
                montageModel->getMontageViewports(getTabNumber(),
                                                  montageViewports);
                
                bool isValid = false;
                bool isFlat  = false;
                bool isLeft = false;
                bool isLateral = true;
                const int32_t numViewports = static_cast<int32_t>(montageViewports.size());
                for (int32_t ivp = 0; ivp < numViewports; ivp++) {
                    const SurfaceMontageViewport& smv = montageViewports[ivp];
                    if (smv.isInside(mousePressX,
                                     mousePressY)) {
                        switch (smv.projectionViewType) {
                            case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_LATERAL:
                                isLeft = true;
                                isLateral = true;
                                break;
                            case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_MEDIAL:
                                isLeft = true;
                                isLateral = false;
                                break;
                            case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_FLAT_SURFACE:
                                isLeft = true;
                                isFlat = true;
                                break;
                            case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_LATERAL:
                                isLeft = false;
                                isLateral = true;
                                break;
                            case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_MEDIAL:
                                isLeft = false;
                                isLateral = false;
                                break;
                            case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_FLAT_SURFACE:
                                isLeft = false;
                                isFlat = true;
                                break;
                        }
                        isValid = true;
                    }
                    
                    if (isValid) {
                        break;
                    }
                }
                
                if (isFlat) {
                    /*
                     * No rotation.
                     */
                    dx = 0.0;
                    dy = 0.0;
                }
                else if (isValid) {
                    if (isLeft == false) {
                        dx = -dx;
                    }
                    if (isLateral == false) {
                        dy = -dy;
                    }
                }
            }
            
            Matrix4x4 rotationMatrix = m_viewingTransformation->getRotationMatrix();
            rotationMatrix.rotateX(-dy);
            rotationMatrix.rotateY(dx);
            m_viewingTransformation->setRotationMatrix(rotationMatrix);
        }
    }
    updateYokedBrowserTabs();
}

/**
 * Apply mouse scaling to the displayed model.
 *
 * @param mouseDX
 *    Change in mouse X coordinate.
 * @param mouseDY
 *    Change in mouse Y coordinate.
 */
void
BrowserTabContent::applyMouseScaling(const int32_t /*mouseDX*/,
                                     const int32_t mouseDY)
{
    if (isVolumeSlicesDisplayed()) {
        float scaling = m_volumeSliceViewingTransformation->getScaling();
        if (mouseDY != 0.0) {
            scaling *= (1.0f + (mouseDY * 0.01));
        }
        if (scaling < 0.01) {
            scaling = 0.01;
        }
        m_volumeSliceViewingTransformation->setScaling(scaling);
    }
    else {
        float scaling = m_viewingTransformation->getScaling();
        if (mouseDY != 0.0) {
            scaling *= (1.0f + (mouseDY * 0.01));
        }
        if (scaling < 0.01) {
            scaling = 0.01;
        }
        m_viewingTransformation->setScaling(scaling);
    }
    updateYokedBrowserTabs();
}

/**
 * Apply mouse translation to the displayed model.
 *
 * @param viewportContent
 *    Content of the viewport.
 * @param mousePressX
 *    X coordinate of where mouse was pressed.
 * @param mousePressY
 *    X coordinate of where mouse was pressed.
 * @param mouseDX
 *    Change in mouse X coordinate.
 * @param mouseDY
 *    Change in mouse Y coordinate.
 */
void
BrowserTabContent::applyMouseTranslation(BrainOpenGLViewportContent* viewportContent,
                                         const int32_t mousePressX,
                                         const int32_t mousePressY,
                                         const int32_t mouseDX,
                                         const int32_t mouseDY)
{
    const int tabIndex = getTabNumber();
    
    if (isVolumeSlicesDisplayed()) {
        const float volumeSliceScaling = m_volumeSliceViewingTransformation->getScaling();
        ModelVolume* modelVolume = getDisplayedVolumeModel();
        VolumeMappableInterface* vf = modelVolume->getUnderlayVolumeFile(tabIndex);
        BoundingBox mybox;
        vf->getVoxelSpaceBoundingBox(mybox);
        float cubesize = std::max(std::max(mybox.getDifferenceX(), mybox.getDifferenceY()), mybox.getDifferenceZ());//factor volume bounding box into slowdown for zoomed in
        float slowdown = 0.005f * cubesize / volumeSliceScaling;//when zoomed in, make the movements slower to match - still changes based on viewport currently
        
        float dx = 0.0;
        float dy = 0.0;
        float dz = 0.0;
        switch (this->getSliceViewPlane())
        {
            case VolumeSliceViewPlaneEnum::ALL:
            {
                int viewport[4];
                viewportContent->getModelViewport(viewport);
                const int32_t halfWidth  = viewport[2] / 2;
                const int32_t halfHeight = viewport[3] / 2;
                const int32_t viewportMousePressedX = mousePressX - viewport[0];
                const int32_t viewportMousePressedY = mousePressY - viewport[1];
                bool isRight  = false;
                bool isTop = false;
                if (viewportMousePressedX > halfWidth) {
                    isRight = true;
                }
                if (viewportMousePressedY > halfHeight) {
                    isTop = true;
                }
                //CaretLogInfo("right: " + AString::fromBool(isRight) + " top: " + AString::fromBool(isTop));
                if (isTop)
                {
                    if (isRight)//coronal
                    {
                        dx = mouseDX * slowdown;
                        dz = mouseDY * slowdown;
                    } else {//parasaggital
                        dy = -mouseDX * slowdown;
                        dz = mouseDY * slowdown;
                    }
                } else {
                    if (isRight)//axial
                    {
                        dx = mouseDX * slowdown;
                        dy = mouseDY * slowdown;
                    }//bottom left has no slice
                }
                break;
            }
            case VolumeSliceViewPlaneEnum::AXIAL:
                dx = mouseDX * slowdown;
                dy = mouseDY * slowdown;
                break;
            case VolumeSliceViewPlaneEnum::CORONAL:
                dx = mouseDX * slowdown;
                dz = mouseDY * slowdown;
                break;
            case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                dy = -mouseDX * slowdown;
                dz = mouseDY * slowdown;
                break;
        }
        
        float translation[3];
        m_volumeSliceViewingTransformation->getTranslation(translation);
        translation[0] += dx;
        translation[1] += dy;
        translation[2] += dz;
        m_volumeSliceViewingTransformation->setTranslation(translation);
    }
    else {
        float dx = mouseDX;
        float dy = mouseDY;
        
        ModelSurfaceMontage* montageModel = getDisplayedSurfaceMontageModel();
        if (montageModel != NULL) {
            std::vector<SurfaceMontageViewport> montageViewports;
            montageModel->getMontageViewports(getTabNumber(),
                                              montageViewports);
            
            bool isValid = false;
            bool isLeft = true;
            bool isLateral = false;
            const int32_t numViewports = static_cast<int32_t>(montageViewports.size());
            for (int32_t ivp = 0; ivp < numViewports; ivp++) {
                const SurfaceMontageViewport& smv = montageViewports[ivp];
                if (smv.isInside(mousePressX,
                                 mousePressY)) {
                    switch (smv.projectionViewType) {
                        case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_LATERAL:
                            isLeft = true;
                            isLateral = true;
                            break;
                        case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_MEDIAL:
                            isLeft = true;
                            isLateral = false;
                            break;
                        case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_FLAT_SURFACE:
                            isLeft = true;
                            isLateral = true;
                            break;
                        case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_LATERAL:
                            isLeft = false;
                            isLateral = true;
                            break;
                        case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_MEDIAL:
                            isLeft = false;
                            isLateral = false;
                            break;
                        case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_FLAT_SURFACE:
                            break;
                            isLeft = false;
                            isLateral = true;
                    }
                    isValid = true;
                }
                
                if (isValid) {
                    break;
                }
            }
            
            if (isValid) {
                if (isLeft == false) {
                    dx = -dx;
                }
                if (isLateral == false) {
                    dx = -dx;
                }
                
                float translation[3];
                m_viewingTransformation->getTranslation(translation);
                translation[0] += dx;
                translation[1] += dy;
                m_viewingTransformation->setTranslation(translation);
            }
        }
        else {
            if (getProjectionViewType() == ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_LATERAL) {
                dx = -dx;
            }
            else if (getProjectionViewType() == ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_FLAT_SURFACE) {
                dx = -dx;
            }
            float translation[3];
            m_viewingTransformation->getTranslation(translation);
            translation[0] += dx;
            translation[1] += dy;
            m_viewingTransformation->setTranslation(translation);
        }
    }
    updateYokedBrowserTabs();
}

/**
 * Get the transformations for drawing a model.
 *
 * @param projectionViewType
 *    Type of projection view
 * @param translationOut
 *    Translation
 * @param rotationMatrixOut
 *    OpenGL rotation matrix.
 * @param scalingOut
 *    Scaling.
 */
void
BrowserTabContent::getTransformationsForOpenGLDrawing(const ProjectionViewTypeEnum::Enum projectionViewType,
                                                      float translationOut[3],
                                                      double rotationMatrixOut[16],
                                                      float& scalingOut) const
{
    if (isVolumeSlicesDisplayed()) {
        m_volumeSliceViewingTransformation->getTranslation(translationOut);
        
        Matrix4x4 rotationMatrix = m_volumeSliceViewingTransformation->getRotationMatrix();
        rotationMatrix.getMatrixForOpenGL(rotationMatrixOut);
        
        scalingOut = m_volumeSliceViewingTransformation->getScaling();
    }
    else {
        m_viewingTransformation->getTranslation(translationOut);

        Matrix4x4 rotationMatrix = m_viewingTransformation->getRotationMatrix();
        double rotationX, rotationY, rotationZ;
        rotationMatrix.getRotation(rotationX,
                                      rotationY,
                                      rotationZ);
        const double rotationFlippedX = -rotationX;
        const double rotationFlippedY = 180.0 - rotationY;
        
        switch (projectionViewType) {
            case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_LATERAL:
                break;
            case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_MEDIAL:
                break;
            case ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_FLAT_SURFACE:
                rotationX =     0.0;
                rotationY =     0.0;
                rotationZ =     0.0;
                break;
            case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_LATERAL:
                rotationX = rotationFlippedX;
                rotationY = rotationFlippedY;
                break;
            case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_MEDIAL:
                rotationX = rotationFlippedX;
                rotationY = rotationFlippedY;
                break;
            case ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_FLAT_SURFACE:
                rotationX =   0.0;
                rotationY = 180.0;
                rotationZ =   0.0;
                break;
        }
        
        Matrix4x4 matrix;
        matrix.setRotation(rotationX,
                           rotationY,
                           rotationZ);
        matrix.getMatrixForOpenGL(rotationMatrixOut);
        
        scalingOut = m_viewingTransformation->getScaling();
    }
}

/**
 * Place the transformations for the given window tab into
 * the model transform.
 * @param windowTabNumber
 *    Tab number for transformations.
 * @param modelTransform
 *    Model transform into which transformations are loaded.
 */
void
BrowserTabContent::getTransformationsInModelTransform(ModelTransform& modelTransform) const
{
    modelTransform.setTranslation(getTranslation());
    
    const Matrix4x4 rotMatrix = getRotationMatrix();
    float m[4][4];
    rotMatrix.getMatrix(m);
    modelTransform.setRotation(m);
    
    const Matrix4x4 obliqueRotationMatrix = getObliqueVolumeRotationMatrix();
    float mob[4][4];
    obliqueRotationMatrix.getMatrix(mob);
    modelTransform.setObliqueRotation(mob);
    
    modelTransform.setScaling(getScaling());
}

/**
 * Apply the transformations to the browser tab.
 * @param modelTransform
 *    Model transform into which transformations are retrieved.
 */
void
BrowserTabContent::setTransformationsFromModelTransform(const ModelTransform& modelTransform)
{
    float translation[3];
    modelTransform.getTranslation(translation);
    const float tx = translation[0];
    const float ty = translation[1];
    const float tz = translation[2];
    
    setTranslation(tx, ty, tz);

    float m[4][4];
    modelTransform.getRotation(m);
    
    Matrix4x4 rotationMatrix;
    rotationMatrix.setMatrix(m);
    setRotationMatrix(rotationMatrix);

    float mob[4][4];
    modelTransform.getObliqueRotation(mob);
    Matrix4x4 obliqueRotationMatrix;
    obliqueRotationMatrix.setMatrix(mob);
    setObliqueVolumeRotationMatrix(obliqueRotationMatrix);

    const float scale = modelTransform.getScaling();
    setScaling(scale);
    updateYokedBrowserTabs();
}


/**
 * Create a scene for an instance of a class.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    saving the scene.
 *
 * @param instanceName
 *    Name of the class' instance.
 *
 * @return Pointer to SceneClass object representing the state of
 *    this object.  Under some circumstances a NULL pointer may be
 *    returned.  Caller will take ownership of returned object.
 */
SceneClass* 
BrowserTabContent::saveToScene(const SceneAttributes* sceneAttributes,
                                const AString& instanceName)
{
    SceneClass* sceneClass = new SceneClass(instanceName,
                                            "BrowserTabContent",
                                            2); // version 2 as of 4/1/2013

    float obliqueMatrix[16];
    m_obliqueVolumeRotationMatrix->getMatrixForOpenGL(obliqueMatrix);
    sceneClass->addFloatArray("m_obliqueVolumeRotationMatrix", obliqueMatrix, 16);
    
    m_sceneClassAssistant->saveMembers(sceneAttributes, 
                                       sceneClass);
    
    return sceneClass;
}

/**
 * Restore the state of an instance of a class.
 * 
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass for the instance of a class that implements
 *     this interface.  May be NULL for some types of scenes.
 */
void 
BrowserTabContent::restoreFromScene(const SceneAttributes* sceneAttributes,
                              const SceneClass* sceneClass)
{
    if (sceneClass == NULL) {
        return;
    }
    
    m_sceneClassAssistant->restoreMembers(sceneAttributes, 
                                          sceneClass);
    
    m_obliqueVolumeRotationMatrix->identity();
    float obliqueMatrix[16];
    const int32_t numInObliqueArray = sceneClass->getFloatArrayValue("m_obliqueVolumeRotationMatrix", obliqueMatrix, 16);
    if (numInObliqueArray == 16) {
        m_obliqueVolumeRotationMatrix->setMatrixFromOpenGL(obliqueMatrix);
    }
    
    /*
     * In older version of workbench, transformation were stored in the
     * model for each tab, so try to restore them.
     */
    if (sceneClass->getVersionNumber() < 2) {
        float translation[3];
        float scaling;
        float rotationMatrix[4][4];
        
        const Model* model = getModelControllerForDisplay();
        if (model != NULL) {
            const bool valid = model->getOldSceneTransformation(m_tabNumber,
                                                                translation,
                                                                scaling,
                                                                rotationMatrix);
            if (valid) {
                setTranslation(translation);
                setScaling(scaling);
                Matrix4x4 m;
                m.setMatrix(rotationMatrix);
                
                ModelSurface* ms = getDisplayedSurfaceModel();
                if (ms != NULL) {
                    /*
                     * Right hemispheres need rotations changed for
                     * proper viewing.
                     */
                    const StructureEnum::Enum structure = ms->getSurface()->getStructure();
                    if (StructureEnum::isRight(structure)) {
                        double rotationX, rotationY, rotationZ;
                        m.getRotation(rotationX,
                                      rotationY,
                                      rotationZ);
                        //rotationX = -rotationX;
                        //rotationY = 180.0 - rotationY;
                        rotationY = 90 + rotationY;
                        rotationZ = -rotationZ;
                        m.identity();
                        m.setRotation(rotationX,
                                           rotationY,
                                           rotationZ);
                    }
                }
                setRotationMatrix(m);
            }
        }
    }
    
    if (getDisplayedWholeBrainModel() != NULL) {
        /*
         * As of 19sep2013 whole brain and volume slice settings were merged
         * (whole brain slice settings removed).  For compatibility, if a
         * whole brain model is being viewed and whole brain slice settings
         * are found, allow them to override volume slice settings.
         */
        const SceneClass* wholeBrainVolumeSettings = sceneClass->getClass("m_wholeBrainSliceSettings");
        if (wholeBrainVolumeSettings != NULL) {
            VolumeSliceSettings settings;
            settings.restoreFromScene(sceneAttributes, wholeBrainVolumeSettings);
            *m_volumeSliceSettings = settings;
        }
    }
}

/**
 * Is the clipping plane enabled?
 * @param indx
 *    Index of clipping plane
 * @return
 *    true if clipping plane enabled, else false.
 */
bool
BrowserTabContent::isClippingPlaneEnabled(const int32_t indx) const
{
    CaretAssertArrayIndex(m_clippingEnabled, 3, indx);
    return m_clippingEnabled[indx];
}

/**
 * Set the clipping plane enabled
 * @param indx
 *    Index of clipping plane
 * status
 *    true if clipping plane enabled, else false.
 */
void
BrowserTabContent::setClippingPlaneEnabled(const int32_t indx,
                                           const bool status)
{
    CaretAssertArrayIndex(m_clippingEnabled, 3, indx);
    m_clippingEnabled[indx] = status;
}

/**
 * Get the clipping plane thickness.
 * @param indx
 *    Index of clipping plane
 * @return
 *    Thickeness of the clipping plane.
 */
float
BrowserTabContent::getClippingPlaneThickness(const int32_t indx) const
{
    CaretAssertArrayIndex(m_clippingThickness, 3, indx);
    return m_clippingThickness[indx];
}

/**
 * Set the clipping plane thickness.
 * @param indx
 *    Index of clipping plane
 * @param value
 *    Thickeness of the clipping plane.
 */
void
BrowserTabContent::setClippingPlaneThickness(const int32_t indx,
                                             const float value)
{
    CaretAssertArrayIndex(m_clippingThickness, 3, indx);
    m_clippingThickness[indx] = value;
}

/**
 * Get the clipping plane coordinate
 * @param indx
 *    Index of clipping plane
 * @return
 *    Coordinate of the clipping plane.
 */
float
BrowserTabContent::getClippingPlaneCoordinate(const int32_t indx) const
{
    CaretAssertArrayIndex(m_clippingCoordinate, 3, indx);
    return m_clippingCoordinate[indx];
}

/**
 * Set the clipping plane coordinate.
 * @param indx
 *    Index of clipping plane
 * @param value
 *    Coordinate of the clipping plane.
 */
void
BrowserTabContent::setClippingPlaneCoordinate(const int32_t indx,
                                const float value)
{
    CaretAssertArrayIndex(m_clippingCoordinate, 3, indx);
    m_clippingCoordinate[indx] = value;
}

/**
 * @return the projection view type (view from left/right)
 */
ProjectionViewTypeEnum::Enum
BrowserTabContent::getProjectionViewType() const
{
    ProjectionViewTypeEnum::Enum projectionViewType = ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_LATERAL;
    
    const ModelSurface* modelSurface = getDisplayedSurfaceModel();
    if (modelSurface != NULL) {
        const SurfaceFile* surfaceFile = modelSurface->getSurface();
        if (surfaceFile != NULL) {
            const StructureEnum::Enum structure = surfaceFile->getStructure();
            const SurfaceTypeEnum::Enum surfaceType = surfaceFile->getSurfaceType();
            if (StructureEnum::isRight(structure)) {
                projectionViewType = ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_LATERAL;
                if (surfaceType == SurfaceTypeEnum::FLAT) {
                    projectionViewType = ProjectionViewTypeEnum::PROJECTION_VIEW_RIGHT_FLAT_SURFACE;
                }
            }
            else {
                projectionViewType = ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_LATERAL;
                if (surfaceType == SurfaceTypeEnum::FLAT) {
                    projectionViewType = ProjectionViewTypeEnum::PROJECTION_VIEW_LEFT_FLAT_SURFACE;
                }
            }
        }
    }
    
    return projectionViewType;
}

/**
 * @return The slice view plane.
 *
 */
VolumeSliceViewPlaneEnum::Enum
BrowserTabContent::getSliceViewPlane() const
{
    return m_volumeSliceSettings->getSliceViewPlane();
}

/**
 * Set the slice view plane.
 * @param windowTabNumber
 *    New value for slice plane.
 */
void
BrowserTabContent::setSliceViewPlane(const VolumeSliceViewPlaneEnum::Enum slicePlane)
{
    m_volumeSliceSettings->setSliceViewPlane(slicePlane);
    updateYokedBrowserTabs();
}

/**
 * @return the slice viewing mode.
 */
VolumeSliceViewModeEnum::Enum
BrowserTabContent::getSliceViewMode() const
{
    return m_volumeSliceSettings->getSliceViewMode();
}

/**
 * Set the slice viewing mode.
 * @param sliceViewMode
 *    New value for view mode
 */
void
BrowserTabContent::setSliceViewMode(const VolumeSliceViewModeEnum::Enum sliceViewMode)
{
    m_volumeSliceSettings->setSliceViewMode(sliceViewMode);
    updateYokedBrowserTabs();
}

/**
 * @return the montage number of columns for the given window tab.
 */
int32_t
BrowserTabContent::getMontageNumberOfColumns() const
{
    return m_volumeSliceSettings->getMontageNumberOfColumns();
}


/**
 * Set the montage number of columns in the given window tab.
 * @param montageNumberOfColumns
 *    New value for montage number of columns
 */
void
BrowserTabContent::setMontageNumberOfColumns(const int32_t montageNumberOfColumns)
{
    m_volumeSliceSettings->setMontageNumberOfColumns(montageNumberOfColumns);
    updateYokedBrowserTabs();
}

/**
 * @return the montage number of rows for the given window tab.
 */
int32_t
BrowserTabContent::getMontageNumberOfRows() const
{
    return m_volumeSliceSettings->getMontageNumberOfRows();
}

/**
 * Set the montage number of rows.
 * @param montageNumberOfRows
 *    New value for montage number of rows
 */
void
BrowserTabContent::setMontageNumberOfRows(const int32_t montageNumberOfRows)
{
    m_volumeSliceSettings->setMontageNumberOfRows(montageNumberOfRows);
    updateYokedBrowserTabs();
}

/**
 * @return the montage slice spacing.
 */
int32_t
BrowserTabContent::getMontageSliceSpacing() const
{
    return m_volumeSliceSettings->getMontageSliceSpacing();
}

/**
 * Set the montage slice spacing.
 * @param montageSliceSpacing
 *    New value for montage slice spacing
 */
void
BrowserTabContent::setMontageSliceSpacing(const int32_t montageSliceSpacing)
{
    m_volumeSliceSettings->setMontageSliceSpacing(montageSliceSpacing);
    updateYokedBrowserTabs();
}

/**
 * Set the selected slices to the origin.
 */
void
BrowserTabContent::setSlicesToOrigin()
{
    selectSlicesAtOrigin();
    updateYokedBrowserTabs();
}

/**
 * Reset the slices.
 */
void
BrowserTabContent::reset()
{
    if (isVolumeSlicesDisplayed()
        || isWholeBrainDisplayed()) {
        m_volumeSliceSettings->reset();
        m_obliqueVolumeRotationMatrix->identity();
    }
    updateYokedBrowserTabs();
}

/**
 * Update the slices coordinates so that they are valid for
 * the given VolumeFile.
 * @param volumeFile
 *   File for which slice coordinates are made valid.
 */
void
BrowserTabContent::updateForVolumeFile(const VolumeMappableInterface* volumeFile)
{
    m_volumeSliceSettings->updateForVolumeFile(volumeFile);
}

/**
 * Set the slice indices so that they are at the origin.
 */
void
BrowserTabContent::selectSlicesAtOrigin()
{
    m_volumeSliceSettings->selectSlicesAtOrigin();
    updateYokedBrowserTabs();
}

/**
 * Set the selected slices to the given coordinate.
 * @param xyz
 *    Coordinate for selected slices.
 */
void
BrowserTabContent::selectSlicesAtCoordinate(const float xyz[3])
{
    m_volumeSliceSettings->selectSlicesAtCoordinate(xyz);
    updateYokedBrowserTabs();
}

/**
 * Return the axial slice index.
 * @return
 *   Axial slice index or negative if invalid
 */
int64_t
BrowserTabContent::getSliceIndexAxial(const VolumeMappableInterface* volumeFile) const
{
    return m_volumeSliceSettings->getSliceIndexAxial(volumeFile);
}

/**
 * Set the axial slice index.
 * @param
 *    New value for axial slice index.
 */
void
BrowserTabContent::setSliceIndexAxial(const VolumeMappableInterface* volumeFile,
                                        const int64_t sliceIndexAxial)
{
    m_volumeSliceSettings->setSliceIndexAxial(volumeFile, sliceIndexAxial);
    updateYokedBrowserTabs();
}

/**
 * Return the coronal slice index.
 * @return
 *   Coronal slice index.
 */
int64_t
BrowserTabContent::getSliceIndexCoronal(const VolumeMappableInterface* volumeFile) const
{
    return m_volumeSliceSettings->getSliceIndexCoronal(volumeFile);
}


/**
 * Set the coronal slice index.
 * @param sliceIndexCoronal
 *    New value for coronal slice index.
 */
void
BrowserTabContent::setSliceIndexCoronal(const VolumeMappableInterface* volumeFile,
                                          const int64_t sliceIndexCoronal)
{
    m_volumeSliceSettings->setSliceIndexCoronal(volumeFile, sliceIndexCoronal);
    updateYokedBrowserTabs();
}

/**
 * Return the parasagittal slice index.
 * @return
 *   Parasagittal slice index.
 */
int64_t
BrowserTabContent::getSliceIndexParasagittal(const VolumeMappableInterface* volumeFile) const
{
    return m_volumeSliceSettings->getSliceIndexParasagittal(volumeFile);
}

/**
 * Set the parasagittal slice index.
 * @param sliceIndexParasagittal
 *    New value for parasagittal slice index.
 */
void
BrowserTabContent::setSliceIndexParasagittal(const VolumeMappableInterface* volumeFile,
                                               const int64_t sliceIndexParasagittal)
{
    m_volumeSliceSettings->setSliceIndexParasagittal(volumeFile,
                                                         sliceIndexParasagittal);
    updateYokedBrowserTabs();
}

/**
 * @return Coordinate of axial slice.
 */
float
BrowserTabContent::getSliceCoordinateAxial() const
{
    return m_volumeSliceSettings->getSliceCoordinateAxial();
}

/**
 * Set the coordinate for the axial slice.
 * @param z
 *    Z-coordinate for axial slice.
 */
void
BrowserTabContent::setSliceCoordinateAxial(const float z)
{
    m_volumeSliceSettings->setSliceCoordinateAxial(z);
    updateYokedBrowserTabs();
}

/**
 * @return Coordinate of coronal slice.
 */
float
BrowserTabContent::getSliceCoordinateCoronal() const
{
    return m_volumeSliceSettings->getSliceCoordinateCoronal();
}

/**
 * Set the coordinate for the coronal slice.
 * @param y
 *    Y-coordinate for coronal slice.
 */
void
BrowserTabContent::setSliceCoordinateCoronal(const float y)
{
    m_volumeSliceSettings->setSliceCoordinateCoronal(y);
    updateYokedBrowserTabs();
}

/**
 * @return Coordinate of parasagittal slice.
 */
float
BrowserTabContent::getSliceCoordinateParasagittal() const
{
    return m_volumeSliceSettings->getSliceCoordinateParasagittal();
}

/**
 * Set the coordinate for the parasagittal slice.
 * @param x
 *    X-coordinate for parasagittal slice.
 */
void
BrowserTabContent::setSliceCoordinateParasagittal(const float x)
{
    m_volumeSliceSettings->setSliceCoordinateParasagittal(x);
    updateYokedBrowserTabs();
}

/**
 * Is the parasagittal slice enabled?
 * @return
 *    Enabled status of parasagittal slice.
 */
bool
BrowserTabContent::isSliceParasagittalEnabled() const
{
    return m_volumeSliceSettings->isSliceParasagittalEnabled();
}

/**
 * Set the enabled status of the parasagittal slice.
 * @param sliceEnabledParasagittal
 *    New enabled status.
 */
void
BrowserTabContent::setSliceParasagittalEnabled(const bool sliceEnabledParasagittal)
{
    m_volumeSliceSettings->setSliceParasagittalEnabled(sliceEnabledParasagittal);
    updateYokedBrowserTabs();
}

/**
 * Is the coronal slice enabled?
 * @return
 *    Enabled status of coronal slice.
 */
bool
BrowserTabContent::isSliceCoronalEnabled() const
{
    return m_volumeSliceSettings->isSliceCoronalEnabled();
}

/**
 * Set the enabled status of the coronal slice.
 * @param sliceEnabledCoronal
 *    New enabled status.
 */
void
BrowserTabContent::setSliceCoronalEnabled(const bool sliceEnabledCoronal)
{
    m_volumeSliceSettings->setSliceCoronalEnabled(sliceEnabledCoronal);
    updateYokedBrowserTabs();
}

/**
 * Is the axial slice enabled?
 * @return
 *    Enabled status of axial slice.
 */
bool
BrowserTabContent::isSliceAxialEnabled() const
{
    return m_volumeSliceSettings->isSliceAxialEnabled();
}

/**
 * Set the enabled status of the axial slice.
 * @param sliceEnabledAxial
 *    New enabled status.
 */
void
BrowserTabContent::setSliceAxialEnabled(const bool sliceEnabledAxial)
{
    m_volumeSliceSettings->setSliceAxialEnabled(sliceEnabledAxial);
    updateYokedBrowserTabs();
}



/**
 * @return Enabled status for left cerebral cortex.
 */
bool
BrowserTabContent::isWholeBrainLeftEnabled() const
{
    return m_wholeBrainSurfaceSettings->isLeftEnabled();
}

/**
 * Set the enabled status for the left hemisphere.
 * @param windowTabNumber
 *    Index of window tab.
 * @param enabled
 *    New enabled status.
 */
void
BrowserTabContent::setWholeBrainLeftEnabled(const bool enabled)
{
    m_wholeBrainSurfaceSettings->setLeftEnabled(enabled);
    updateYokedBrowserTabs();
}

/**
 * @return Enabled status for right cerebral cortex.
 */
bool
BrowserTabContent::isWholeBrainRightEnabled() const
{
    return m_wholeBrainSurfaceSettings->isRightEnabled();
}

/**
 * Set the enabled status for the right hemisphere.
 * @param enabled
 *    New enabled status.
 */
void
BrowserTabContent::setWholeBrainRightEnabled(const bool enabled)
{
    m_wholeBrainSurfaceSettings->setRightEnabled(enabled);
    updateYokedBrowserTabs();
}

/**
 * @return Enabled status for cerebellum.
 */
bool
BrowserTabContent::isWholeBrainCerebellumEnabled() const
{
    return m_wholeBrainSurfaceSettings->isCerebellumEnabled();
}

/**
 * Set the enabled status for the cerebellum.
 * @param enabled
 *    New enabled status.
 */
void
BrowserTabContent::setWholeBrainCerebellumEnabled(const bool enabled)
{
    m_wholeBrainSurfaceSettings->setCerebellumEnabled(enabled);
    updateYokedBrowserTabs();
}

/**
 * @return The separation between the left and right surfaces.
 */
float
BrowserTabContent::getWholeBrainLeftRightSeparation() const
{
    return m_wholeBrainSurfaceSettings->getLeftRightSeparation();
}

/**
 * Set the separation between the cerebellum and the left/right surfaces.
 * @param separation
 *     New value for separation.
 */
void
BrowserTabContent::setWholeBrainLeftRightSeparation(const float separation)
{
    m_wholeBrainSurfaceSettings->setLeftRightSeparation(separation);
    updateYokedBrowserTabs();
}

/**
 * @return The separation between the left/right surfaces.
 */
float
BrowserTabContent::getWholeBrainCerebellumSeparation() const
{
    return m_wholeBrainSurfaceSettings->getCerebellumSeparation();
}

/**
 * Set the separation between the cerebellum and the eft and right surfaces.
 * @param separation
 *     New value for separation.
 */
void
BrowserTabContent::setWholeBrainCerebellumSeparation(const float separation)
{
    m_wholeBrainSurfaceSettings->setCerebellumSeparation(separation);
    updateYokedBrowserTabs();
}

/**
 * @return Selected yoking group.
 */
YokingGroupEnum::Enum
BrowserTabContent::getYokingGroup() const
{
    return m_yokingGroup;
}

/**
 * Set the selected yoking group.
 *
 * @param yokingGroup
 *    New value for yoking group.
 */
void
BrowserTabContent::setYokingGroup(const YokingGroupEnum::Enum yokingGroup)
{
    m_yokingGroup = yokingGroup;
    
    if (m_yokingGroup == YokingGroupEnum::YOKING_GROUP_OFF) {
        return;
    }
    
    /*
     * Find another browser tab using the same yoking as 'me' and copy
     * yoked data from the other browser tab.
     */
    for (std::set<BrowserTabContent*>::iterator iter = s_allBrowserTabContent.begin();
         iter != s_allBrowserTabContent.end();
         iter++) {
        BrowserTabContent* btc = *iter;
        if (btc != this) {
            if (btc->getYokingGroup() == m_yokingGroup) {
                *m_viewingTransformation = *btc->m_viewingTransformation;
                *m_volumeSliceViewingTransformation = *btc->m_volumeSliceViewingTransformation;
                *m_volumeSliceSettings = *btc->m_volumeSliceSettings;
                break;
            }
        }
    }
}

/**
 * @return Is this browser tab yoked?
 */
bool
BrowserTabContent::isYoked() const
{
    const bool yoked = (m_yokingGroup != YokingGroupEnum::YOKING_GROUP_OFF);
    return yoked;
}

/**
 * Update other browser tabs with yoked data.
 */
void
BrowserTabContent::updateYokedBrowserTabs()
{
    if (m_yokingGroup == YokingGroupEnum::YOKING_GROUP_OFF) {
        return;
    }
    
    /*
     * Copy yoked data from 'me' to all other yoked browser tabs
     */
    for (std::set<BrowserTabContent*>::iterator iter = s_allBrowserTabContent.begin();
         iter != s_allBrowserTabContent.end();
         iter++) {
        BrowserTabContent* btc = *iter;
        if (btc != this) {
            if (btc->getYokingGroup() == m_yokingGroup) {
                *btc->m_viewingTransformation = *m_viewingTransformation;
                *btc->m_volumeSliceViewingTransformation = *m_volumeSliceViewingTransformation;
                *btc->m_volumeSliceSettings = *m_volumeSliceSettings;
                *btc->m_obliqueVolumeRotationMatrix = *m_obliqueVolumeRotationMatrix;
                //*btc->m_wholeBrainSurfaceSettings = *m_wholeBrainSurfaceSettings;
            }
        }
    }
}


