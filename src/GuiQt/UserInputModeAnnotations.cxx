
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

#define __USER_INPUT_MODE_ANNOTATIONS_DECLARE__
#include "UserInputModeAnnotations.h"
#undef __USER_INPUT_MODE_ANNOTATIONS_DECLARE__

#include <cmath>

#include "AnnotationBrowserTab.h"
#include "AnnotationChangeCoordinateDialog.h"
#include "AnnotationCreateDialog.h"
#include "AnnotationColorBar.h"
#include "AnnotationCoordinate.h"
#include "AnnotationCoordinateInformation.h"
#include "AnnotationFile.h"
#include "AnnotationManager.h"
#include "AnnotationMultiCoordinateShape.h"
#include "AnnotationTwoCoordinateShape.h"
#include "AnnotationPasteDialog.h"
#include "AnnotationRedoUndoCommand.h"
#include "AnnotationSpatialModification.h"
#include "AnnotationText.h"
#include "AnnotationTextEditorDialog.h"
#include "AnnotationOneCoordinateShape.h"
#include "Brain.h"
#include "BrainBrowserWindow.h"
#include "BrainOpenGLViewportContent.h"
#include "BrainOpenGLWidget.h"
#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "CaretUndoStack.h"
#include "CursorEnum.h"
#include "CaretPreferences.h"
#include "DisplayPropertiesAnnotation.h"
#include "EventAnnotationCreateNewType.h"
#include "EventAnnotationGetDrawnInWindow.h"
#include "EventGraphicsUpdateAllWindows.h"
#include "EventIdentificationRequest.h"
#include "EventUserInterfaceUpdate.h"
#include "EventManager.h"
#include "EventGraphicsUpdateOneWindow.h"
#include "GestureEvent.h"
#include "GuiManager.h"
#include "IdentificationManager.h"
#include "KeyEvent.h"
#include "ModelSurfaceMontage.h"
#include "MouseEvent.h"
#include "SelectionItemAnnotation.h"
#include "SelectionManager.h"
#include "SelectionItemSurfaceNode.h"
#include "SelectionItemVoxel.h"
#include "SessionManager.h"
#include "Surface.h"
#include "SurfaceMontageViewport.h"
#include "UserInputModeAnnotationsContextMenu.h"
#include "UserInputModeAnnotationsWidget.h"
#include "WuQMessageBox.h"

using namespace caret;


/**
 * \class caret::UserInputModeAnnotations 
 * \brief Input mode processor for Annotations.
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param windowIndex
 *     Index of window
 */
UserInputModeAnnotations::UserInputModeAnnotations(const int32_t windowIndex)
: UserInputModeAnnotations(UserInputModeEnum::Enum::ANNOTATIONS,
                           windowIndex)
{
}

/**
 * Constructor.
 *
 * @param userInputMode
 *     Input
 * @param windowIndex
 *     Index of window
 */
UserInputModeAnnotations::UserInputModeAnnotations(const UserInputModeEnum::Enum userInputMode,
                                                   const int32_t windowIndex)
: UserInputModeView(windowIndex,
                    userInputMode),
m_browserWindowIndex(windowIndex),
m_annotationUnderMouse(NULL),
m_annotationBeingDragged(NULL)
{
    m_allowMultipleSelectionModeFlag = true;
    m_mode = MODE_SELECT;
    m_annotationUnderMouseSizeHandleType = AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE;
    m_annotationUnderMousePolyLineCoordinateIndex = -1;
    
    m_modeNewAnnotationFileSpaceAndType.grabNew(new NewAnnotationFileSpaceAndType(NULL,
                                                                                  AnnotationCoordinateSpaceEnum::VIEWPORT,
                                                                                  AnnotationTypeEnum::LINE));
    m_newAnnotationCreatingWithMouseDrag.grabNew(NULL);
    
    m_annotationToolsWidget = new UserInputModeAnnotationsWidget(this,
                                                                 m_browserWindowIndex);
    setWidgetForToolBar(m_annotationToolsWidget);
    
    EventManager::get()->addEventListener(this, EventTypeEnum::EVENT_ANNOTATION_CREATE_NEW_TYPE);
}

/**
 * Destructor.
 */
UserInputModeAnnotations::~UserInputModeAnnotations()
{
    EventManager::get()->removeAllEventsFromListener(this);
}

/**
 * Receive an event.
 *
 * @param event
 *    An event for which this instance is listening.
 */
void
UserInputModeAnnotations::receiveEvent(Event* event)
{
    if (event->getEventType() == EventTypeEnum::EVENT_ANNOTATION_CREATE_NEW_TYPE) {
        EventAnnotationCreateNewType* annotationEvent = dynamic_cast<EventAnnotationCreateNewType*>(event);
        CaretAssert(annotationEvent);
        
        AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
        
        annotationManager->deselectAllAnnotationsForEditing(m_browserWindowIndex);
        resetAnnotationUnderMouse();
        
        const AnnotationTypeEnum::Enum annType(annotationEvent->getAnnotationType());
        m_modeNewAnnotationFileSpaceAndType.grabNew(new NewAnnotationFileSpaceAndType(annotationEvent->getAnnotationFile(),
                                                                                      annotationEvent->getAnnotationSpace(),
                                                                                      annType));

        Mode mode(MODE_NEW_WITH_DRAG_START);
        switch (annType) {
            case AnnotationTypeEnum::BOX:
                break;
            case AnnotationTypeEnum::BROWSER_TAB:
                break;
            case AnnotationTypeEnum::COLOR_BAR:
                break;
            case AnnotationTypeEnum::IMAGE:
                break;
            case AnnotationTypeEnum::LINE:
                break;
            case AnnotationTypeEnum::OVAL:
                break;
            case AnnotationTypeEnum::POLY_LINE:
                switch (annotationEvent->getPolyLineDrawingMode()) {
                    case EventAnnotationCreateNewType::CONTINUOUS:
                        mode = MODE_NEW_WITH_DRAG_START;
                        break;
                    case EventAnnotationCreateNewType::DISCRETE:
                        mode = MODE_NEW_WITH_CLICK_SERIES_START;
                        break;
                }
                break;
            case AnnotationTypeEnum::SCALE_BAR:
                break;
            case AnnotationTypeEnum::TEXT:
                break;
        }
        
        setMode(mode);
        EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    }
}


/**
 * Called when 'this' user input receiver is set
 * to receive events.
 */
void
UserInputModeAnnotations::initialize()
{
    m_mode = MODE_SELECT;
    DisplayPropertiesAnnotation* dpa = GuiManager::get()->getBrain()->getDisplayPropertiesAnnotation();
    dpa->setDisplayAnnotations(true);
    resetAnnotationUnderMouse();
}

/**
 * Called when 'this' user input receiver is no
 * longer set to receive events.
 */
void
UserInputModeAnnotations::finish()
{
    m_mode = MODE_SELECT;
    resetAnnotationUnderMouse();
}

void
UserInputModeAnnotations::resetAnnotationUnderMouse()
{
    m_annotationUnderMouse  = NULL;
    m_annotationUnderMouseSizeHandleType = AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE;
    m_annotationUnderMousePolyLineCoordinateIndex = -1;
    m_annotationBeingDragged = NULL;
    m_annotationBeingDraggedHandleType = AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE;
}

/**
 * Called to update the input receiver for various events.
 */
void
UserInputModeAnnotations::update()
{
    
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString
UserInputModeAnnotations::toString() const
{
    return "UserInputModeBorders";
}

/**
 * @return the mode.
 */
UserInputModeAnnotations::Mode
UserInputModeAnnotations::getMode() const
{
    return m_mode;
}

/**
 * Set the mode.
 * @param mode
 *    New value for mode.
 */
void
UserInputModeAnnotations::setMode(const Mode mode)
{
    if (m_mode != mode) {
        m_mode = mode;
        
        bool drawingModeFlag(false);
        switch (m_mode) {
            case MODE_NEW_WITH_DRAG_START:
                drawingModeFlag = true;
                break;
            case MODE_NEW_WITH_CLICK_SERIES:
                drawingModeFlag = true;
                break;
            case MODE_NEW_WITH_CLICK_SERIES_START:
                drawingModeFlag = true;
                break;
            case MODE_NEW_WITH_DRAG:
                drawingModeFlag = true;
                break;
            case MODE_PASTE:
                break;
            case MODE_PASTE_SPECIAL:
                break;
            case MODE_SELECT:
                break;
            case MODE_SET_COORDINATE_ONE:
                break;
            case MODE_SET_COORDINATE_TWO:
                break;
        }
        
        if ( ! drawingModeFlag) {
            resetAnnotationBeingCreated();
        }
    }
    EventManager::get()->sendEvent(EventUserInterfaceUpdate().getPointer());
}

/**
 * @return The cursor for display in the OpenGL widget.
 */
CursorEnum::Enum
UserInputModeAnnotations::getCursor() const
{
    
    CursorEnum::Enum cursor = CursorEnum::CURSOR_DEFAULT;
    
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_NEW_WITH_DRAG:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_PASTE:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_PASTE_SPECIAL:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_SELECT:
            if (m_annotationUnderMouse != NULL) {
                cursor = CursorEnum::CURSOR_FOUR_ARROWS;
                
                switch (m_annotationUnderMouseSizeHandleType) {
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM:
                        cursor = CursorEnum::CURSOR_RESIZE_VERTICAL;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_LEFT:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_LEFT_TOP_RIGHT;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_BOTTOM_RIGHT:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_RIGHT_TOP_LEFT;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_LEFT:
                        cursor = CursorEnum::CURSOR_RESIZE_HORIZONTAL;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_RIGHT:
                        cursor = CursorEnum::CURSOR_RESIZE_HORIZONTAL;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP:
                        cursor = CursorEnum::CURSOR_RESIZE_VERTICAL;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_LEFT:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_RIGHT_TOP_LEFT;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_BOX_TOP_RIGHT:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_LEFT_TOP_RIGHT;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_END:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_LEFT_TOP_RIGHT;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_LINE_START:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_LEFT_TOP_RIGHT;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE:
                        cursor = CursorEnum::CURSOR_FOUR_ARROWS;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION:
                        cursor = CursorEnum::CURSOR_ROTATION;
                        break;
                    case AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_POLY_LINE_COORDINATE:
                        cursor = CursorEnum::CURSOR_RESIZE_BOTTOM_LEFT_TOP_RIGHT;
                        break;
                }
            }
            break;
        case MODE_SET_COORDINATE_ONE:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
        case MODE_SET_COORDINATE_TWO:
            cursor = CursorEnum::CURSOR_CROSS;
            break;
    }
    
    return cursor;
}

/**
 * Delete all selected annotations except color bars which are turned off for display
 */
void
UserInputModeAnnotations::deleteSelectedAnnotations()
{
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    if (annotationManager->isAnnotationSelectedForEditingDeletable(m_browserWindowIndex)) {
        std::vector<Annotation*> selectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(m_browserWindowIndex);
        if ( ! selectedAnnotations.empty()) {
            AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
            undoCommand->setModeDeleteAnnotations(selectedAnnotations);
            AString errorMessage;
            if ( !  annotationManager->applyCommand(getUserInputMode(),
                                                    undoCommand,
                                                    errorMessage)) {
                WuQMessageBox::errorOk(m_annotationToolsWidget,
                                       errorMessage);
            }
        }
    }
    else {
        std::vector<Annotation*> selectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(m_browserWindowIndex);
        for (std::vector<Annotation*>::iterator iter = selectedAnnotations.begin();
             iter != selectedAnnotations.end();
             iter++) {
            Annotation* ann = *iter;
            CaretAssert(ann);
            
            if (ann->getType() == AnnotationTypeEnum::COLOR_BAR) {
                AnnotationColorBar* colorBar = dynamic_cast<AnnotationColorBar*>(ann);
                CaretAssert(colorBar);
                colorBar->setDisplayed(false);
            }
        }
    }

    EventManager::get()->sendEvent(EventUserInterfaceUpdate().getPointer());
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
}

/**
 * Process a key press event
 *
 * @param keyEvent
 *     Key event information.
 * @return
 *     True if the input process recognized the key event
 *     and the key event SHOULD NOT be propagated to parent
 *     widgets
 */
bool
UserInputModeAnnotations::keyPressEvent(const KeyEvent& keyEvent)
{
    bool keyWasProcessedFlag(false);
    
    const int32_t keyCode = keyEvent.getKeyCode();
    switch (keyCode) {
        case Qt::Key_Backspace:
        case Qt::Key_Delete:
        {
            switch (m_mode) {
                case MODE_NEW_WITH_DRAG_START:
                    break;
                case MODE_NEW_WITH_CLICK_SERIES:
                    break;
                case MODE_NEW_WITH_CLICK_SERIES_START:
                    break;
                case MODE_NEW_WITH_DRAG:
                    break;
                case MODE_PASTE:
                    break;
                case MODE_PASTE_SPECIAL:
                    break;
                case MODE_SELECT:
                    deleteSelectedAnnotations();
                    keyWasProcessedFlag = true;
                    break;
                case MODE_SET_COORDINATE_ONE:
                    break;
                case MODE_SET_COORDINATE_TWO:
                    break;
            }
        }
            break;
        case Qt::Key_Escape:
        {
            bool selectModeFlag = false;
            switch (m_mode) {
                case MODE_NEW_WITH_DRAG_START:
                    break;
                case MODE_NEW_WITH_CLICK_SERIES:
                    selectModeFlag = true;
                    break;
                case MODE_NEW_WITH_CLICK_SERIES_START:
                    selectModeFlag = true;
                    break;
                case MODE_NEW_WITH_DRAG:
                    break;
                case MODE_PASTE:
                    selectModeFlag = true;
                    break;
                case MODE_PASTE_SPECIAL:
                    selectModeFlag = true;
                    break;
                case MODE_SELECT:
                    break;
                case MODE_SET_COORDINATE_ONE:
                    selectModeFlag = true;
                    break;
                case MODE_SET_COORDINATE_TWO:
                    selectModeFlag = true;
                    break;
            }
            if (selectModeFlag) {
                setMode(MODE_SELECT);
                EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
                keyWasProcessedFlag = true;
            }
        }
            break;
        case Qt::Key_Down:
        case Qt::Key_Left:
        case Qt::Key_Right:
        case Qt::Key_Up:
        {
            Annotation* selectedAnnotation = getSingleSelectedAnnotation();
            if (selectedAnnotation != NULL) {
                bool changeCoordFlag  = false;
                bool moveOnePixelFlag = false;
                switch (selectedAnnotation->getCoordinateSpace()) {
                    case AnnotationCoordinateSpaceEnum::CHART:
                        changeCoordFlag = true;
                        break;
                    case AnnotationCoordinateSpaceEnum::SPACER:
                        changeCoordFlag = true;
                        moveOnePixelFlag = true;
                        break;
                    case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
                        changeCoordFlag = true;
                        break;
                    case AnnotationCoordinateSpaceEnum::SURFACE:
                        break;
                    case AnnotationCoordinateSpaceEnum::TAB:
                        changeCoordFlag  = true;
                        moveOnePixelFlag = true;
                        break;
                    case AnnotationCoordinateSpaceEnum::VIEWPORT:
                        break;
                    case AnnotationCoordinateSpaceEnum::WINDOW:
                        changeCoordFlag  = true;
                        moveOnePixelFlag = true;
                        break;
                }
                
                if (changeCoordFlag) {
                    keyWasProcessedFlag = true;
                    
                    float distanceX = 1.0;
                    float distanceY = 1.0;
                    if (moveOnePixelFlag) {
                        const float pixelHeight = keyEvent.getOpenGLWidget()->height();
                        const float pixelWidth  = keyEvent.getOpenGLWidget()->width();
                        /*
                         * 100 is "full width/height" in relative coordinates.
                         */
                        distanceX = 100.0 / pixelWidth;
                        distanceY = 100.0 / pixelHeight;
                    }
                    if (keyEvent.isShiftKeyDownFlag()) {
                        const float multiplier = 10;
                        distanceX *= multiplier;
                        distanceY *= multiplier;
                    }
                    
                    float dx = 0.0;
                    float dy = 0.0;
                    switch (keyCode) {
                        case Qt::Key_Down:
                            dy = -distanceY;
                            break;
                        case Qt::Key_Left:
                            dx = -distanceX;
                            break;
                        case Qt::Key_Right:
                            dx = distanceX;
                            break;
                        case Qt::Key_Up:
                            dy = distanceY;
                            break;
                        default:
                            CaretAssert(0);
                            break;
                    }
                    
                    AnnotationTwoCoordinateShape* twoCoordShape = selectedAnnotation->castToTwoCoordinateShape();
                    AnnotationOneCoordinateShape* oneCoordShape = selectedAnnotation->castToOneCoordinateShape();
                    AnnotationMultiCoordinateShape* multiCoordShape = selectedAnnotation->castToMultiCoordinateShape();
                    
                    {
                            bool surfaceFlag = false;
                            switch (selectedAnnotation->getCoordinateSpace()) {
                                case AnnotationCoordinateSpaceEnum::CHART:
                                    break;
                                case AnnotationCoordinateSpaceEnum::SPACER:
                                    break;
                                case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
                                    break;
                                case AnnotationCoordinateSpaceEnum::SURFACE:
                                    surfaceFlag = true;
                                    break;
                                case AnnotationCoordinateSpaceEnum::TAB:
                                    break;
                                case AnnotationCoordinateSpaceEnum::VIEWPORT:
                                    break;
                                case AnnotationCoordinateSpaceEnum::WINDOW:
                                    break;
                            }
                            
                            if ( ! surfaceFlag) {
                                AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
                                std::vector<Annotation*> annotations;
                                annotations.push_back(selectedAnnotation);
                                
                                if (twoCoordShape != NULL) {
                                    AnnotationCoordinate startCoord = *twoCoordShape->getStartCoordinate();
                                    float xyzStart[3];
                                    startCoord.getXYZ(xyzStart);
                                    xyzStart[0] += dx;
                                    xyzStart[1] += dy;
                                    startCoord.setXYZ(xyzStart);
                                    
                                    AnnotationCoordinate endCoord   = *twoCoordShape->getEndCoordinate();
                                    float xyzEnd[3];
                                    endCoord.getXYZ(xyzEnd);
                                    xyzEnd[0] += dx;
                                    xyzEnd[1] += dy;
                                    endCoord.setXYZ(xyzEnd);
                                    
                                    undoCommand->setModeCoordinateOneAndTwo(startCoord, endCoord, annotations);
                                }
                                else if (oneCoordShape != NULL) {
                                    AnnotationCoordinate coord = *oneCoordShape->getCoordinate();
                                    float xyz[3];
                                    coord.getXYZ(xyz);
                                    xyz[0] += dx;
                                    xyz[1] += dy;
                                    coord.setXYZ(xyz);
                                    
                                    undoCommand->setModeCoordinateOne(coord,
                                                                      annotations);
                                }
                                else if (multiCoordShape != NULL) {
                                    std::vector<std::unique_ptr<AnnotationCoordinate>> allCoords;
                                    std::vector<std::unique_ptr<const AnnotationCoordinate>> constCoords;
                                    multiCoordShape->getCopyOfAllCoordinates(allCoords);
                                    for (const auto& ac : allCoords) {
                                        float xyz[3];
                                        ac->getXYZ(xyz);
                                        xyz[0] += dx;
                                        xyz[1] += dy;
                                        ac->setXYZ(xyz);
                                        std::unique_ptr<const AnnotationCoordinate> acCopy(new AnnotationCoordinate(*ac));
                                        constCoords.push_back(std::move(acCopy));
                                    }
                                    
                                    undoCommand->setModeCoordinateMulti(constCoords,
                                                                        annotations);
                                }
                                else {
                                    CaretAssert(0);
                                }
                                
                                if ( ! keyEvent.isFirstKeyPressFlag()) {
                                    undoCommand->setMergeEnabled(true);
                                }
                                
                                AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
                                AString errorMessage;
                                if ( ! annMan->applyCommand(getUserInputMode(),
                                                            undoCommand,
                                                            errorMessage)) {
                                    WuQMessageBox::errorOk(m_annotationToolsWidget,
                                                           errorMessage);
                                }
                                
                                EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
                                EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
                            }
                            
                        }
                }
            }
        }
            break;
    }
    
    return keyWasProcessedFlag;
}

/**
 * Initialize user drawing a new annotation.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::initializeUserDrawingNewAnnotation(const MouseEvent& mouseEvent)
{
    if (m_newAnnotationCreatingWithMouseDrag != NULL) {
        m_newAnnotationCreatingWithMouseDrag.grabNew(NULL);
    }
    
    /*
     * Note ALWAYS use WINDOW space for the drag anntotion.
     * Otherwise it will not get displayed if surface/stereotaxic
     */
    m_newAnnotationCreatingWithMouseDrag.grabNew(new NewMouseDragCreateAnnotation(m_modeNewAnnotationFileSpaceAndType->m_annotationFile,
                                                                                  AnnotationCoordinateSpaceEnum::WINDOW,
                                                                                  m_modeNewAnnotationFileSpaceAndType->m_annotationType,
                                                                                  mouseEvent));
}

/**
 * Process a mouse left drag with no keys down event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftDrag(const MouseEvent& mouseEvent)
{
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
        {
            initializeUserDrawingNewAnnotation(mouseEvent);
            m_mode = MODE_NEW_WITH_DRAG;
            return;
        }
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            return;
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            return;
            break;
        case MODE_NEW_WITH_DRAG:
            userDrawingAnnotationFromMouseDrag(mouseEvent);
            return;
            break;
        case MODE_PASTE:
            break;
        case MODE_PASTE_SPECIAL:
            break;
        case MODE_SELECT:
            break;
        case MODE_SET_COORDINATE_ONE:
            break;
        case MODE_SET_COORDINATE_TWO:
            break;
    }
    
    AnnotationCoordinateSpaceEnum::Enum draggingCoordinateSpace = AnnotationCoordinateSpaceEnum::VIEWPORT;
    
    std::vector<Annotation*> selectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(m_browserWindowIndex);
    const int32_t numSelectedAnnotations = static_cast<int32_t>(selectedAnnotations.size());
    
    bool draggingValid = false;
    if (numSelectedAnnotations == 1) {
        draggingCoordinateSpace = selectedAnnotations[0]->getCoordinateSpace();
        draggingValid = true;
    }
    else if (numSelectedAnnotations > 1) {
        if (m_annotationBeingDraggedHandleType == AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE) {
            bool allSameSpaceFlag = true;
            draggingCoordinateSpace = selectedAnnotations[0]->getCoordinateSpace();
            for (int32_t i = 1; i < numSelectedAnnotations; i++) {
                if (selectedAnnotations[i]->getCoordinateSpace() != draggingCoordinateSpace) {
                    allSameSpaceFlag = false;
                    break;
                }
            }
            
            if (allSameSpaceFlag) {
                draggingValid = true;
            }
        }
    }
    
    switch (draggingCoordinateSpace) {
        case AnnotationCoordinateSpaceEnum::CHART:
            break;
        case AnnotationCoordinateSpaceEnum::SPACER:
            break;
        case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            break;
        case AnnotationCoordinateSpaceEnum::SURFACE:
            break;
        case AnnotationCoordinateSpaceEnum::TAB:
            break;
        case AnnotationCoordinateSpaceEnum::VIEWPORT:
            /*
             * No dragging in viewport space
             */
            draggingValid = false;
            break;
        case AnnotationCoordinateSpaceEnum::WINDOW:
            break;
    }
    
    if (draggingValid) {
        BrainOpenGLViewportContent* vpContent = mouseEvent.getViewportContent();
        if (vpContent == NULL) {
            return;
        }
        
        
        AnnotationCoordinateInformation coordInfo;
        AnnotationCoordinateInformation::createCoordinateInformationFromXY(mouseEvent.getOpenGLWidget(),
                                                                 mouseEvent.getViewportContent(),
                                                                 mouseEvent.getX(),
                                                                 mouseEvent.getY(),
                                                                 coordInfo);
        
        float spaceOriginX = 0.0;
        float spaceOriginY = 0.0;
        float spaceWidth  = 0.0;
        float spaceHeight = 0.0;
        
        switch (draggingCoordinateSpace) {
            case AnnotationCoordinateSpaceEnum::CHART:
            {
                /*
                 * Want viewport within montage that contains the surface
                 */
                Matrix4x4 projectionMatrix;
                Matrix4x4 modelviewMatrix;
                int32_t chartVP[4];
                if (vpContent->getChartDataMatricesAndViewport(projectionMatrix,
                                                               modelviewMatrix,
                                                               chartVP)) {
                    spaceOriginX = chartVP[0];
                    spaceOriginY = chartVP[1];
                    spaceWidth   = chartVP[2];
                    spaceHeight  = chartVP[3];
                }
            }
                break;
            case AnnotationCoordinateSpaceEnum::SPACER:
            {
                int viewport[4];
                vpContent->getTabViewportBeforeApplyingMargins(viewport);
                spaceOriginX = viewport[0];
                spaceOriginY = viewport[1];
                spaceWidth   = viewport[2];
                spaceHeight  = viewport[3];
            }
                break;
            case AnnotationCoordinateSpaceEnum::STEREOTAXIC:
            {
                /*
                 * Want viewport within montage that contains the surface
                 */
                int montVP[4];
                vpContent->getSurfaceMontageModelViewport(mouseEvent.getX(),
                                                          mouseEvent.getY(),
                                                          montVP);
                spaceOriginX = montVP[0];
                spaceOriginY = montVP[1];
                spaceWidth   = montVP[2];
                spaceHeight  = montVP[3];
            }
                break;
            case AnnotationCoordinateSpaceEnum::SURFACE:
            {
                /*
                 * Want viewport within montage that contains the surface
                 */
                int montVP[4];
                vpContent->getSurfaceMontageModelViewport(mouseEvent.getX(),
                                                          mouseEvent.getY(),
                                                          montVP);
                spaceOriginX = montVP[0];
                spaceOriginY = montVP[1];
                spaceWidth   = montVP[2];
                spaceHeight  = montVP[3];
                
                /*
                 * Also need to find viewport where mouse was originally
                 * pressed.  If the viewports are different, then 
                 * mouse has left the viewport in which it was pressed
                 * and the mouse position is no longer valid.
                 */
                int pressVP[4];
                vpContent->getSurfaceMontageModelViewport(mouseEvent.getPressedX(),
                                                          mouseEvent.getPressedY(),
                                                          pressVP);
                if ((montVP[0] != pressVP[0])
                    || (montVP[1] != pressVP[1])
                    || (montVP[2] != pressVP[2])
                    || (montVP[3] != pressVP[3])) {
                    draggingValid = false;
                }
            }
                break;
            case AnnotationCoordinateSpaceEnum::TAB:
            {
                int viewport[4];
                vpContent->getModelViewport(viewport);
                spaceOriginX = viewport[0];
                spaceOriginY = viewport[1];
                spaceWidth   = viewport[2];
                spaceHeight  = viewport[3];
            }
                break;
            case AnnotationCoordinateSpaceEnum::VIEWPORT:
                CaretAssert(0);
                break;
            case AnnotationCoordinateSpaceEnum::WINDOW:
            {
                int viewport[4];
                vpContent->getWindowViewport(viewport);
                spaceOriginX = viewport[0];
                spaceOriginY = viewport[1];
                spaceWidth   = viewport[2];
                spaceHeight  = viewport[3];
            }
                break;
        }
        
        if (draggingValid) {
            const float dx = mouseEvent.getDx();
            const float dy = mouseEvent.getDy();
            
            const float mouseViewportX = mouseEvent.getX() - spaceOriginX;
            const float mouseViewportY = mouseEvent.getY() - spaceOriginY;
            
            const float mousePressViewportX = mouseEvent.getPressedX() - spaceOriginX;
            const float mousePressViewportY = mouseEvent.getPressedY() - spaceOriginY;
            
            AnnotationSpatialModification annSpatialMod(m_annotationBeingDraggedHandleType,
                                                        spaceWidth,
                                                        spaceHeight,
                                                        mousePressViewportX,
                                                        mousePressViewportY,
                                                        mouseViewportX,
                                                        mouseViewportY,
                                                        dx,
                                                        dy,
                                                        m_annotationUnderMousePolyLineCoordinateIndex,
                                                        mouseEvent.isFirstDragging());
            if (coordInfo.m_surfaceSpaceInfo.m_validFlag) {
                annSpatialMod.setSurfaceCoordinateAtMouseXY(coordInfo.m_surfaceSpaceInfo.m_structure,
                                                            coordInfo.m_surfaceSpaceInfo.m_numberOfNodes,
                                                            coordInfo.m_surfaceSpaceInfo.m_nodeIndex);
            }
            
            if (coordInfo.m_modelSpaceInfo.m_validFlag) {
                annSpatialMod.setStereotaxicCoordinateAtMouseXY(coordInfo.m_modelSpaceInfo.m_xyz[0],
                                                                coordInfo.m_modelSpaceInfo.m_xyz[1],
                                                                coordInfo.m_modelSpaceInfo.m_xyz[2]);
            }
            
            if (coordInfo.m_chartSpaceInfo.m_validFlag) {
                annSpatialMod.setChartCoordinateAtMouseXY(coordInfo.m_chartSpaceInfo.m_xyz[0],
                                                          coordInfo.m_chartSpaceInfo.m_xyz[1],
                                                          coordInfo.m_chartSpaceInfo.m_xyz[2]);
            }
            
            if ((dx != 0.0)
                || (dy != 0.0)) {
                AnnotationCoordinateInformation previousMouseXYCoordInfo;
                AnnotationCoordinateInformation::createCoordinateInformationFromXY(mouseEvent.getOpenGLWidget(),
                                                                                   mouseEvent.getViewportContent(),
                                                                                   mouseEvent.getX() - dx,
                                                                                   mouseEvent.getY() - dy,
                                                                                   previousMouseXYCoordInfo);
                if (previousMouseXYCoordInfo.m_chartSpaceInfo.m_validFlag) {
                    annSpatialMod.setChartCoordinateAtPreviousMouseXY(previousMouseXYCoordInfo.m_chartSpaceInfo.m_xyz[0],
                                                                      previousMouseXYCoordInfo.m_chartSpaceInfo.m_xyz[1],
                                                                      previousMouseXYCoordInfo.m_chartSpaceInfo.m_xyz[2]);
                }
            }
            
            std::vector<Annotation*> annotationsBeforeMoveAndResize;
            std::vector<Annotation*> annotationsAfterMoveAndResize;
            
            for (int32_t i = 0; i < numSelectedAnnotations; i++) {
                Annotation* annotationModified(selectedAnnotations[i]->clone());
                if (annotationModified->applySpatialModification(annSpatialMod)) {
                    annotationsBeforeMoveAndResize.push_back(selectedAnnotations[i]);
                    annotationsAfterMoveAndResize.push_back(annotationModified);
                }
                else {
                    delete annotationModified;
                    annotationModified = NULL;
                }
            }
            CaretAssert(annotationsAfterMoveAndResize.size() == annotationsBeforeMoveAndResize.size());
            
            if ( ! annotationsAfterMoveAndResize.empty()) {
                AnnotationRedoUndoCommand* command = new AnnotationRedoUndoCommand();
                command->setModeLocationAndSize(annotationsBeforeMoveAndResize,
                                                annotationsAfterMoveAndResize);
                
                if ( ! mouseEvent.isFirstDragging()) {
                    command->setMergeEnabled(true);
                }
                
                AString errorMessage;
                if ( !  annotationManager->applyCommand(getUserInputMode(),
                                                        command,
                                                        errorMessage)) {
                    WuQMessageBox::errorOk(m_annotationToolsWidget,
                                           errorMessage);
                }
            }
            
            for (std::vector<Annotation*>::iterator iter = annotationsAfterMoveAndResize.begin();
                 iter != annotationsAfterMoveAndResize.end();
                 iter++) {
                delete *iter;
            }
            annotationsAfterMoveAndResize.clear();
            
            EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
            EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
        }
    }
}

/**
 * Process a mouse left drag with only the alt key down event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftDragWithAlt(const MouseEvent& /*mouseEvent*/)
{
}

/**
 * Process a mouse left drag with ctrl key down event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftDragWithCtrl(const MouseEvent& /*mouseEvent*/) { }

/**
 * Process a mouse left drag with ctrl and shift keys down event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftDragWithCtrlShift(const MouseEvent& /*mouseEvent*/) { }

/**
 * Process a mouse left drag with shift key down event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftDragWithShift(const MouseEvent& mouseEvent)
{
    mouseLeftDrag(mouseEvent);
}

/**
 * Process a mouse left click event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftClick(const MouseEvent& mouseEvent)
{
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            processModeNewMouseLeftClick(mouseEvent);
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            userDrawingAnnotationFromMouseDrag(mouseEvent);
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
        {
            initializeUserDrawingNewAnnotation(mouseEvent);
            m_mode = MODE_NEW_WITH_CLICK_SERIES;
            AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
            annotationManager->setAnnotationBeingDrawnInWindow(m_browserWindowIndex,
                                                               m_newAnnotationCreatingWithMouseDrag->getAnnotation());
            EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
        }
            break;
        case MODE_NEW_WITH_DRAG:
            break;
        case MODE_PASTE:
            pasteAnnotationFromAnnotationClipboard(mouseEvent);
            break;
        case MODE_PASTE_SPECIAL:
            pasteAnnotationFromAnnotationClipboardAndChangeSpace(mouseEvent);
            break;
        case MODE_SELECT:
            break;
        case MODE_SET_COORDINATE_ONE:
            processModeSetCoordinate(mouseEvent);
            break;
        case MODE_SET_COORDINATE_TWO:
            processModeSetCoordinate(mouseEvent);
            break;
    }
}

/**
 * Process a mouse left click with Shift key event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftClickWithShift(const MouseEvent& mouseEvent)
{
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            /*
             * Insert vertex and finish annotation
             */
            userDrawingAnnotationFromMouseDrag(mouseEvent);
            createNewAnnotationFromMouseDrag(mouseEvent);
            m_mode = MODE_SELECT;
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            WuQMessageBox::errorOk(m_annotationToolsWidget,
                                   "Annotation has not been started.  "
                                   "Click mouse WITHOUT SHIFT key down to draw annotation or press ESC key to exit drawing.");
            break;
        case MODE_NEW_WITH_DRAG:
            break;
        case MODE_PASTE:
            break;
        case MODE_PASTE_SPECIAL:
            break;
        case MODE_SELECT:
            if (m_allowMultipleSelectionModeFlag) {
                const bool shiftKeyDown(true);
                const bool singleSelectionModeFlag(false);
                processMouseSelectAnnotation(mouseEvent,
                                             shiftKeyDown,
                                             singleSelectionModeFlag);
            }
            break;
        case MODE_SET_COORDINATE_ONE:
            break;
        case MODE_SET_COORDINATE_TWO:
            break;
    }
}

/**
 * Process a mouse left press event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftPress(const MouseEvent& mouseEvent)
{
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            break;
        case MODE_NEW_WITH_DRAG:
            break;
        case MODE_PASTE:
            break;
        case MODE_PASTE_SPECIAL:
            break;
        case MODE_SELECT:
        {
            /*
             * Single click selects clicked annotation and deselects any that are selected
             */
            const bool shiftKeyDown(false);
            const bool singleSelectionModeFlag(true);
            processMouseSelectAnnotation(mouseEvent,
                                         shiftKeyDown,
                                         singleSelectionModeFlag);
        }
            break;
        case MODE_SET_COORDINATE_ONE:
            break;
        case MODE_SET_COORDINATE_TWO:
            break;
    }
}

/**
 * Set the annotation under the mouse that results in update of the cursor.
 *
 * @param mouseEvent
 *     Mouse event information.
 * @param annotationIDIn
 *     Optional (if not NULL) annotation ID that may be provided by caller
 *     and can be used to avoid an identification operation.
 */
void
UserInputModeAnnotations::setAnnotationUnderMouse(const MouseEvent& mouseEvent,
                                                  SelectionItemAnnotation* annotationIDIn)
{
    m_annotationUnderMouseSizeHandleType = AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_NONE;
    
    
    BrainOpenGLWidget* openGLWidget = mouseEvent.getOpenGLWidget();
    SelectionItemAnnotation* annotationID = annotationIDIn;
    if (annotationID == NULL) {
        annotationID = openGLWidget->performIdentificationAnnotations(mouseEvent.getX(),
                                                                      mouseEvent.getY());
    }

    if (annotationID->isValid()) {
        m_annotationUnderMouse = annotationID->getAnnotation();
        m_annotationUnderMouseSizeHandleType = annotationID->getSizingHandle();
        m_annotationUnderMousePolyLineCoordinateIndex = annotationID->getPolyLineCoordinateIndex();
    }
    else {
        m_annotationUnderMouse = NULL; 
    }
    
    openGLWidget->updateCursor();
    
    /*
     * Update graphics only if an annotation was passed to this method (WB-820)
     */
    if (annotationIDIn != NULL) {
#if BRAIN_OPENGL_INFO_SUPPORTS_DISPLAY_LISTS
        openGLWidget->updateGL();
#else
        openGLWidget->update();
#endif
    }
}

/**
 * User drawing a new annotation from dragging the mouse from corner/end to another
 * corner/end.
 *
 * @param mouseEvent
 *     Mouse event issued when mouse button was released.
 */
void
UserInputModeAnnotations::userDrawingAnnotationFromMouseDrag(const MouseEvent& mouseEvent)
{
    if (m_newAnnotationCreatingWithMouseDrag != NULL) {
        m_newAnnotationCreatingWithMouseDrag->update(mouseEvent,
                                                     mouseEvent.getX(),
                                                     mouseEvent.getY());

        AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
        annotationManager->setAnnotationBeingDrawnInWindow(m_browserWindowIndex,
                                                           m_newAnnotationCreatingWithMouseDrag->getAnnotation());
        EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    }
}

/**
 * Create a new annotation from dragging the mouse from corner/end to another
 * corner/end.
 *
 * @param mouseEvent
 *     Mouse event issued when mouse button was released.
 */
void
UserInputModeAnnotations::createNewAnnotationFromMouseDrag(const MouseEvent& mouseEvent)
{
    if (m_newAnnotationCreatingWithMouseDrag != NULL) {
        
        switch (m_mode) {
            case MODE_NEW_WITH_CLICK_SERIES:
                break;
            case MODE_NEW_WITH_CLICK_SERIES_START:
                break;
            case MODE_NEW_WITH_DRAG:
                break;
            case MODE_NEW_WITH_DRAG_START:
                break;
            case MODE_PASTE:
                break;
            case MODE_PASTE_SPECIAL:
                break;
            case MODE_SELECT:
                break;
            case MODE_SET_COORDINATE_ONE:
                break;
            case MODE_SET_COORDINATE_TWO:
                break;
        }
        std::vector<Vector3D> coords = m_newAnnotationCreatingWithMouseDrag->getDrawingCoordinates();
        
        Annotation* ann = AnnotationCreateDialog::newAnnotationFromSpaceTypeAndBounds(mouseEvent,
                                                                                      coords,
                                                                                      m_modeNewAnnotationFileSpaceAndType->m_annotationSpace,
                                                                                      m_modeNewAnnotationFileSpaceAndType->m_annotationType,
                                                                                      m_modeNewAnnotationFileSpaceAndType->m_annotationFile);
        if (ann != NULL) {
            selectAnnotation(ann);
        }
        
        setMode(MODE_SELECT);
        
        EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    }
}

/**
 * Reset the annotation that is being created
 */
void
UserInputModeAnnotations::resetAnnotationBeingCreated()
{
    m_newAnnotationCreatingWithMouseDrag.grabNew(NULL);
    
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    annotationManager->setAnnotationBeingDrawnInWindow(m_browserWindowIndex,
                                                       NULL);
}

/**
 * Process a mouse left release event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftRelease(const MouseEvent& mouseEvent)
{
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            break;
        case MODE_NEW_WITH_DRAG:
            createNewAnnotationFromMouseDrag(mouseEvent);
            m_mode = MODE_SELECT;
            break;
        case MODE_PASTE:
            break;
        case MODE_PASTE_SPECIAL:
            break;
        case MODE_SELECT:
            break;
        case MODE_SET_COORDINATE_ONE:
            break;
        case MODE_SET_COORDINATE_TWO:
            break;
    }
    
    m_annotationBeingDragged = NULL;
    
    setAnnotationUnderMouse(mouseEvent,
                            NULL);
}

/**
 * Process a mouse left double-click event.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseLeftDoubleClick(const MouseEvent& mouseEvent)
{
    const int32_t mouseX = mouseEvent.getX();
    const int32_t mouseY = mouseEvent.getY();
    
    BrainOpenGLWidget* openGLWidget = mouseEvent.getOpenGLWidget();
    SelectionItemAnnotation* annotationID = openGLWidget->performIdentificationAnnotations(mouseX,
                                                                            mouseY);
    if (annotationID->isValid()) {
        Annotation* annotation = annotationID->getAnnotation();
        if (annotation != NULL) {
            AnnotationText* textAnnotation = dynamic_cast<AnnotationText*>(annotation);
            if (textAnnotation != NULL) {
                AnnotationTextEditorDialog ted(textAnnotation,
                                               openGLWidget);
                /*
                 * Note: Y==0 is at top for widget.
                 *       Y==0 is at bottom for OpenGL mouse x,y
                 */
                ted.move(openGLWidget->mapToGlobal(QPoint(mouseX,
                                                          (openGLWidget->height() - mouseY + 20))));
                ted.exec();
                EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
            }
        }
    }
}


/**
 * Process a mouse move with no buttons or keys down
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseMove(const MouseEvent& mouseEvent)
{
    setAnnotationUnderMouse(mouseEvent,
                            NULL);
}

/**
 * Process a mouse move with no buttons and shift key down
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::mouseMoveWithShift(const MouseEvent& mouseEvent)
{
    setAnnotationUnderMouse(mouseEvent,
                            NULL);
}

/**
 * Process a gesture event (pinch zoom; or rotate)
 *
 * @param gestureEvent
 *     Gesture event information.
 */
void
UserInputModeAnnotations::gestureEvent(const GestureEvent& gestureEvent)
{
    BrainOpenGLViewportContent* viewportContent = gestureEvent.getViewportContent();
    if (viewportContent == NULL) {
        return;
    }
    
    BrowserTabContent* browserTabContent = viewportContent->getBrowserTabContent();
    if (browserTabContent == NULL) {
        return;
    }
    
    switch (gestureEvent.getType()) {
        case GestureEvent::Type::PINCH:
            break;
        case GestureEvent::Type::ROTATE:
        {
            float deltaRotateAngle = gestureEvent.getValue();
            
            AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
            std::vector<Annotation*> selectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(m_browserWindowIndex);
            
            float rotationAngle(0.0);
            std::vector<AnnotationOneCoordinateShape*> twoDimAnns;
            for (auto a : selectedAnnotations) {
                AnnotationOneCoordinateShape* a2d = a->castToOneCoordinateShape();
                if (a2d != NULL) {
                    if (a2d->isSizeHandleValid(AnnotationSizingHandleTypeEnum::ANNOTATION_SIZING_HANDLE_ROTATION)) {
                        if (twoDimAnns.empty()) {
                            rotationAngle = a2d->getRotationAngle() + deltaRotateAngle;
                        }
                        twoDimAnns.push_back(a2d);
                    }
                }
            }
            
            if ( ! twoDimAnns.empty()) {
                AnnotationRedoUndoCommand* command = new AnnotationRedoUndoCommand();
                std::vector<Annotation*> modAnns(twoDimAnns.begin(),
                                                 twoDimAnns.end());
                command->setModeRotationAngle(rotationAngle,
                                              modAnns);
                AString errorMessage;
                if ( !  annotationManager->applyCommand(getUserInputMode(),
                                                        command,
                                                        errorMessage)) {
                    WuQMessageBox::errorOk(m_annotationToolsWidget,
                                           errorMessage);
                }
                EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
                EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
            }
        }
    }
}

/**
 * Process a mouse left click to set a coordinate.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::processModeSetCoordinate(const MouseEvent& mouseEvent)
{
    Annotation* selectedAnnotation = getSingleSelectedAnnotation();
    if (selectedAnnotation == NULL) {
        return;
    }

    AnnotationCoordinateInformation coordInfo;
    AnnotationCoordinateInformation::createCoordinateInformationFromXY(mouseEvent.getOpenGLWidget(),
                                                             mouseEvent.getViewportContent(),
                                                             mouseEvent.getX(),
                                                             mouseEvent.getY(),
                                                             coordInfo);
    
    AnnotationTwoCoordinateShape* twoCoordShape = selectedAnnotation->castToTwoCoordinateShape();
    AnnotationOneCoordinateShape* oneCoordShape = selectedAnnotation->castToOneCoordinateShape();

    AnnotationCoordinate* coordinate = NULL;
    AnnotationCoordinate* otherCoordinate = NULL;
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            break;
        case MODE_NEW_WITH_DRAG:
            break;
        case MODE_PASTE:
            break;
        case MODE_PASTE_SPECIAL:
            break;
        case MODE_SELECT:
            break;
        case MODE_SET_COORDINATE_ONE:
            if (twoCoordShape != NULL) {
                coordinate      = twoCoordShape->getStartCoordinate();
                otherCoordinate = twoCoordShape->getEndCoordinate();
            }
            else if (oneCoordShape != NULL) {
                coordinate = oneCoordShape->getCoordinate();
            }
            break;
        case MODE_SET_COORDINATE_TWO:
            if (twoCoordShape != NULL) {
                coordinate      = twoCoordShape->getEndCoordinate();
                otherCoordinate = twoCoordShape->getStartCoordinate();
            }
            break;
    }
    
    if (coordinate != NULL) {
        StructureEnum::Enum structure = StructureEnum::INVALID;
        int32_t numNodes = -1;
        int32_t nodeIndex = -1;
        float surfaceOffset = 0.0;
        AnnotationSurfaceOffsetVectorTypeEnum::Enum surfaceVectorType = AnnotationSurfaceOffsetVectorTypeEnum::CENTROID_THRU_VERTEX;
        coordinate->getSurfaceSpace(structure, numNodes, nodeIndex, surfaceOffset, surfaceVectorType);
        coordInfo.m_surfaceSpaceInfo.m_nodeOffsetLength = surfaceOffset;
        coordInfo.m_surfaceSpaceInfo.m_nodeVectorOffsetType = surfaceVectorType;
    }
    
    AnnotationChangeCoordinateDialog changeCoordDialog(coordInfo,
                                                       selectedAnnotation,
                                                       coordinate,
                                                       otherCoordinate,
                                                       m_annotationToolsWidget);
    if (changeCoordDialog.exec() == AnnotationChangeCoordinateDialog::Accepted) {
        
    }

    setMode(MODE_SELECT);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
}


/**
 * Process a mouse left click for new mode.
 *
 * @param mouseEvent
 *     Mouse event information.
 */
void
UserInputModeAnnotations::processModeNewMouseLeftClick(const MouseEvent& mouseEvent)
{
    resetAnnotationUnderMouse();
    
    std::vector<Vector3D> coords;
    coords.emplace_back(mouseEvent.getPressedX(),
                        mouseEvent.getPressedY(),
                        0.0);
    Annotation* ann = AnnotationCreateDialog::newAnnotationFromSpaceAndType(mouseEvent,
                                                                            coords,
                                                                            m_modeNewAnnotationFileSpaceAndType->m_annotationSpace,
                                                                            m_modeNewAnnotationFileSpaceAndType->m_annotationType,
                                                                            m_modeNewAnnotationFileSpaceAndType->m_annotationFile);
    if (ann != NULL) {
        selectAnnotation(ann);
    }
    
    setMode(MODE_SELECT);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
}

/**
 * @return The single selected annotation.  If there is ONE and ONLY one annotation
 * selected, it is returned.  Otherwise, NULL is returned.
 */
Annotation*
UserInputModeAnnotations::getSingleSelectedAnnotation()
{
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    std::vector<Annotation*> allSelectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(m_browserWindowIndex);
    Annotation* selectedAnnotation = NULL;
    if (allSelectedAnnotations.size() == 1) {
        CaretAssertVectorIndex(allSelectedAnnotations, 0);
        selectedAnnotation = allSelectedAnnotations[0];
    }
    return selectedAnnotation;
}


/**
 * Select the given annotation (typically when a new annontation is created).
 *
 * @param annotation
 *    Annotation that is selected.
 */
void
UserInputModeAnnotations::selectAnnotation(Annotation* annotation)
{
    resetAnnotationUnderMouse();
    
    m_annotationBeingDragged = annotation;
    m_annotationUnderMouse   = annotation;
}


/**
 * Process a mouse left click for selection mode.
 *
 * @param mouseEvent
 *     Mouse event information.
 * @param shiftKeyDownFlag
 *     True if shift key is down.
 * @param singleSelectionModeFlag
 *     If true, deselect any other annotations so that only the annotation under mouse is selected
 */
void
UserInputModeAnnotations::processMouseSelectAnnotation(const MouseEvent& mouseEvent,
                                                       const bool shiftKeyDownFlag,
                                                       const bool singleSelectionModeFlag)
{
    BrainOpenGLWidget* openGLWidget = mouseEvent.getOpenGLWidget();
    const int mouseX = mouseEvent.getX();
    const int mouseY = mouseEvent.getY();
    
    m_annotationBeingDragged = NULL;
    
    /*
     * NOTE: When selecting annotations:
     *    (A) When the mouse is clicked WITHOUT the SHIFT key down, the user is in
     *        'single-annotation-selection-mode' and at most one annotation will
     *        be selected when this method completes.
     *    (B) If the mouse is clicked with the SHIFT key down, the user is in
     *        'multi-annotation-selection-mode' and any number of annotation will
     *        be selected when this method completes.
     */
    SelectionItemAnnotation* annotationID = openGLWidget->performIdentificationAnnotations(mouseX,
                                                                                           mouseY);
    Annotation* selectedAnnotation = annotationID->getAnnotation();
    
    /*
     * If only one annotation may be selected, deselect all other annotations
     */
    if (singleSelectionModeFlag) {
        GuiManager::get()->getBrain()->getAnnotationManager()->deselectAllAnnotationsForEditing(m_browserWindowIndex);
    }
    
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    AnnotationManager::SelectionMode selectionMode = AnnotationManager::SELECTION_MODE_SINGLE;
    if (m_allowMultipleSelectionModeFlag) {
        selectionMode = AnnotationManager::SELECTION_MODE_EXTENDED;
    }
    annotationManager->selectAnnotationForEditing(m_browserWindowIndex,
                                        selectionMode,
                                        shiftKeyDownFlag,
                                        selectedAnnotation);
    
    setAnnotationUnderMouse(mouseEvent,
                            annotationID);

    if (selectedAnnotation != NULL) {
        m_annotationBeingDragged = selectedAnnotation;
        m_annotationBeingDraggedHandleType = annotationID->getSizingHandle();
        m_annotationUnderMousePolyLineCoordinateIndex = annotationID->getPolyLineCoordinateIndex();
    }
    
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
}

/**
 * Show a context menu (pop-up menu at mouse location)
 *
 * @param mouseEvent
 *     Mouse event information.
 * @param menuPosition
 *     Point at which menu is displayed (passed to QMenu::exec())
 * @param openGLWidget
 *     OpenGL widget in which context menu is requested
 */
void
UserInputModeAnnotations::showContextMenu(const MouseEvent& mouseEvent,
                                          const QPoint& menuPosition,
                                          BrainOpenGLWidget* openGLWidget)
{
    BrainOpenGLViewportContent* viewportContent = mouseEvent.getViewportContent();
    BrowserTabContent* tabContent = viewportContent->getBrowserTabContent();
    if (tabContent == NULL) {
        return;
    }

    /*
     * Select any annotation that is under the mouse.
     * There might not be an annotation under the
     * mouse and that is okay.
     */
    const bool shiftKeyDown(false);
    const bool singleSelectionModeFlag(false);
    processMouseSelectAnnotation(mouseEvent,
                                 shiftKeyDown,
                                 singleSelectionModeFlag);
    
    UserInputModeAnnotationsContextMenu contextMenu(this,
                                                    mouseEvent,
                                                    NULL,
                                                    tabContent,
                                                    openGLWidget);
    if (contextMenu.actions().size() > 0) {
        contextMenu.exec(menuPosition);
        
        Annotation* newAnnotation = contextMenu.getNewAnnotationCreatedByContextMenu();
        if (newAnnotation != NULL) {
            selectAnnotation(newAnnotation);

            EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
            EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
        }
    }
}

/**
 * @return True if the edit menu items should be enabled, else false.
 */
bool
UserInputModeAnnotations::isEditMenuValid() const
{
    bool editMenuValid = false;
    
    switch (m_mode) {
        case MODE_NEW_WITH_DRAG_START:
            break;
        case MODE_NEW_WITH_CLICK_SERIES:
            break;
        case MODE_NEW_WITH_CLICK_SERIES_START:
            break;
        case MODE_NEW_WITH_DRAG:
            break;
        case MODE_PASTE:
            editMenuValid = true;
            break;
        case MODE_PASTE_SPECIAL:
            editMenuValid = true;
            break;
        case MODE_SELECT:
            editMenuValid = true;
            break;
        case MODE_SET_COORDINATE_ONE:
            break;
        case MODE_SET_COORDINATE_TWO:
            break;
    }
    
    return editMenuValid;
}

/**
 * Cut the selected annotation.  Only functions if
 * one and only one annotation is selected.
 */
void
UserInputModeAnnotations::cutAnnotation()
{
    std::vector<std::pair<Annotation*, AnnotationFile*> > selectedAnnotations;
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    annotationManager->getAnnotationsAndFilesSelectedForEditing(m_browserWindowIndex,
                                                                selectedAnnotations);
    
    if (selectedAnnotations.size() == 1) {
        CaretAssertVectorIndex(selectedAnnotations, 0);
        AnnotationFile* annotationFile = selectedAnnotations[0].second;
        Annotation* annotation = selectedAnnotations[0].first;
        
        annotationManager->copyAnnotationToClipboard(annotationFile,
                                                     annotation);
        
        std::vector<Annotation*> annotationVector;
        annotationVector.push_back(annotation);
        AnnotationRedoUndoCommand* undoCommand = new AnnotationRedoUndoCommand();
        undoCommand->setModeCutAnnotations(annotationVector);
        AString errorMessage;
        if ( ! annotationManager->applyCommand(getUserInputMode(),
                                               undoCommand,
                                               errorMessage)) {
            WuQMessageBox::errorOk(m_annotationToolsWidget,
                                   errorMessage);
        }
        
        EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
        EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    }
}


/**
 * Process a selection that was made from the browser window's edit menu.
 * Intended for override by sub-classes.
 *
 * @param editMenuItem
 *     Item that was selected from the edit menu.
 */
void
UserInputModeAnnotations::processEditMenuItemSelection(const BrainBrowserWindowEditMenuItemEnum::Enum editMenuItem)
{
    if ( ! isEditMenuValid()) {
        return;
    }
    
    
    switch (editMenuItem) {
        case BrainBrowserWindowEditMenuItemEnum::COPY:
        {
            AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
            std::vector<std::pair<Annotation*, AnnotationFile*> > selectedAnnotations;
            annotationManager->getAnnotationsAndFilesSelectedForEditing(m_browserWindowIndex,
                                                      selectedAnnotations);
            
            if (selectedAnnotations.size() == 1) {
                CaretAssertVectorIndex(selectedAnnotations, 0);
                annotationManager->copyAnnotationToClipboard(selectedAnnotations[0].second,
                                                             selectedAnnotations[0].first);
            }
        }
            break;
        case BrainBrowserWindowEditMenuItemEnum::CUT:
            cutAnnotation();
            break;
        case BrainBrowserWindowEditMenuItemEnum::DELETER:
            deleteSelectedAnnotations();
            break;
        case BrainBrowserWindowEditMenuItemEnum::DESELECT_ALL:
            processDeselectAllAnnotations();
            break;
        case BrainBrowserWindowEditMenuItemEnum::PASTE:
        {
            const MouseEvent* mouseEvent = getMousePosition();
            if (mouseEvent != NULL) {
                pasteAnnotationFromAnnotationClipboard(*mouseEvent);
            }
            else {
                setMode(MODE_PASTE);
            }
        }
            break;
        case BrainBrowserWindowEditMenuItemEnum::PASTE_SPECIAL:
        {
            const MouseEvent* mouseEvent = getMousePosition();
            if (mouseEvent != NULL) {
                pasteAnnotationFromAnnotationClipboardAndChangeSpace(*mouseEvent);
            }
            else {
                setMode(MODE_PASTE_SPECIAL);
            }
        }
            break;
        case BrainBrowserWindowEditMenuItemEnum::REDO:
        {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
            CaretUndoStack* undoStack = annMan->getCommandRedoUndoStack(getUserInputMode());
            
            AString errorMessage;
            if ( ! undoStack->redo(errorMessage)) {
                WuQMessageBox::errorOk(m_annotationToolsWidget,
                                       errorMessage);
            }
            
            EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
            EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
        }
            break;
        case BrainBrowserWindowEditMenuItemEnum::SELECT_ALL:
            processSelectAllAnnotations();
            break;
        case BrainBrowserWindowEditMenuItemEnum::UNDO:
        {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
            CaretUndoStack* undoStack = annMan->getCommandRedoUndoStack(getUserInputMode());
            
            AString errorMessage;
            if ( ! undoStack->undo(errorMessage)) {
                WuQMessageBox::errorOk(m_annotationToolsWidget,
                                       errorMessage);
            }
            
            EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
            EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
        }
            break;
    }
}

/**
 * Process the deselection of all annotations.
 */
void
UserInputModeAnnotations::processDeselectAllAnnotations()
{
    std::vector<Annotation*> annotationsSelected;
    
    switch (getUserInputMode()) {
        case UserInputModeEnum::Enum::ANNOTATIONS:
        {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
            annMan->deselectAllAnnotationsForEditing(m_browserWindowIndex);
        }
            break;
        case UserInputModeEnum::Enum::BORDERS:
        case UserInputModeEnum::Enum::FOCI:
        case UserInputModeEnum::Enum::IMAGE:
        case UserInputModeEnum::Enum::INVALID:
            break;
        case UserInputModeEnum::Enum::TILE_TABS_MANUAL_LAYOUT_EDITING:
        {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
            annMan->deselectAllAnnotationsForEditing(m_browserWindowIndex);
        }
            break;
        case UserInputModeEnum::Enum::VIEW:
            break;
        case UserInputModeEnum::Enum::VOLUME_EDIT:
            break;
    }
        
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
}

/**
 * Process the selection of all annotations.
 */
void
UserInputModeAnnotations::processSelectAllAnnotations()
{
    std::vector<Annotation*> annotationsSelected;
    
    switch (getUserInputMode()) {
        case UserInputModeEnum::Enum::ANNOTATIONS:
        {
            EventAnnotationGetDrawnInWindow getDrawnEvent(m_browserWindowIndex);
            EventManager::get()->sendEvent(getDrawnEvent.getPointer());
            getDrawnEvent.getAnnotations(annotationsSelected);

            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
            annMan->deselectAllAnnotationsForEditing(m_browserWindowIndex);

            annMan->setAnnotationsForEditing(m_browserWindowIndex,
                                             annotationsSelected);
        }
            break;
        case UserInputModeEnum::Enum::BORDERS:
        case UserInputModeEnum::Enum::FOCI:
        case UserInputModeEnum::Enum::IMAGE:
        case UserInputModeEnum::Enum::INVALID:
            break;
        case UserInputModeEnum::Enum::TILE_TABS_MANUAL_LAYOUT_EDITING:
        {
            AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
            annMan->deselectAllAnnotationsForEditing(m_browserWindowIndex);

            BrainBrowserWindow* bbw = GuiManager::get()->getBrowserWindowByWindowIndex(m_browserWindowIndex);
            CaretAssert(bbw);
            std::vector<BrowserTabContent*> allTabContent;
            bbw->getAllTabContent(allTabContent);
            
            std::vector<Annotation*> annotations;
            for (auto btc : allTabContent) {
                AnnotationBrowserTab* bta = btc->getManualLayoutBrowserTabAnnotation();
                if (bta->isBrowserTabDisplayed()) {
                    annotations.push_back(bta);
                }
            }
            
            if ( ! annotations.empty()) {
                AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
                CaretAssert(annMan);
                annMan->setAnnotationsForEditing(m_browserWindowIndex,
                                                 annotations);
            }
            EventManager::get()->sendEvent(EventGraphicsUpdateOneWindow(m_browserWindowIndex).getPointer());
        }
            break;
        case UserInputModeEnum::Enum::VIEW:
            break;
        case UserInputModeEnum::Enum::VOLUME_EDIT:
            break;
    }
    
    
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
}


/**
 * Get the menu items that should be enabled for the current user input processor.
 * Intended for override by sub-classes.
 * Unless this method is overridden, all items on Edit menu are disabled.
 *
 * @param enabledEditMenuItemsOut
 *     Upon exit contains edit menu items that should be enabled.
 * @param redoMenuItemSuffixTextOut
 *     If the redo menu is enabled, the contents of string becomes
 *     the suffix for the 'Redo' menu item.
 * @param undoMenuItemSuffixTextOut
 *     If the undo menu is enabled, the contents of string becomes
 *     the suffix for the 'Undo' menu item.
 * @param pasteTextOut
 *     If not empty, this text is shown for the PASTE menu item
 * @param pasteSpecialTextOut
 *     If not empty, this text is shown for the PASTE_SPECIAL menu item
 */
void
UserInputModeAnnotations::getEnabledEditMenuItems(std::vector<BrainBrowserWindowEditMenuItemEnum::Enum>& enabledEditMenuItemsOut,
                                                  AString& redoMenuItemSuffixTextOut,
                                                  AString& undoMenuItemSuffixTextOut,
                                                  AString& pasteTextOut,
                                                  AString& pasteSpecialTextOut)
{
    enabledEditMenuItemsOut.clear();
    redoMenuItemSuffixTextOut = "";
    undoMenuItemSuffixTextOut = "";
    
    pasteTextOut        = BrainBrowserWindowEditMenuItemEnum::toGuiName(BrainBrowserWindowEditMenuItemEnum::PASTE);
    pasteSpecialTextOut = BrainBrowserWindowEditMenuItemEnum::toGuiName(BrainBrowserWindowEditMenuItemEnum::PASTE_SPECIAL);
    
    
    if (isEditMenuValid()) {
        AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
        const std::vector<Annotation*> allSelectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(m_browserWindowIndex);
        
        /*
         * Copy, Cut, and Delete disabled if ANY colorbar is selected.
         */
        bool allAllowCopyCutPasteFlag = true;
        bool allAllowDeleteFlag       = true;
        bool allAllowSelectFlag       = true;
        for (std::vector<Annotation*>::const_iterator iter = allSelectedAnnotations.begin();
             iter != allSelectedAnnotations.end();
             iter++) {
            const Annotation* ann = *iter;
            if ( ! ann->testProperty(Annotation::Property::COPY_CUT_PASTE)) {
                allAllowCopyCutPasteFlag = false;
            }
            
            if (ann->getType() != AnnotationTypeEnum::BROWSER_TAB) {
                if ( ! ann->testProperty(Annotation::Property::DELETION)) {
                    allAllowDeleteFlag = false;
                }
            }
        }
        
        const bool anySelectedFlag = ( ! allSelectedAnnotations.empty());
        const bool oneSelectedFlag = (allSelectedAnnotations.size() == 1);
        if (allAllowCopyCutPasteFlag) {
            if (oneSelectedFlag) {
                enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::COPY);
                enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::CUT);
            }
        }
        
        if (allAllowDeleteFlag) {
            if (anySelectedFlag) {
                enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::DELETER);
            }
        }
        
        if (allAllowSelectFlag) {
            enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::SELECT_ALL);
        }
        if (anySelectedFlag) {
            enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::DESELECT_ALL);
        }
        
        if (annotationManager->isAnnotationOnClipboardValid()) {
            const Annotation* clipBoardAnn = annotationManager->getAnnotationOnClipboard();
            CaretAssert(clipBoardAnn);

            enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::PASTE);
            enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::PASTE_SPECIAL);
            
            clipBoardAnn->getTextForPasteMenuItems(pasteTextOut,
                                                   pasteSpecialTextOut);
        }
        
        CaretUndoStack* undoStack = annotationManager->getCommandRedoUndoStack(getUserInputMode());
        
        if (undoStack->canRedo()) {
            enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::REDO);
            redoMenuItemSuffixTextOut = undoStack->redoText();
        }
        
        if (undoStack->canUndo()) {
            enabledEditMenuItemsOut.push_back(BrainBrowserWindowEditMenuItemEnum::UNDO);
            undoMenuItemSuffixTextOut = undoStack->undoText();
        }
    }
}

/**
 * Paste the annotation from the annotation clipboard.
 *
 * @param mouseEvent
 *     Mouse event containing XY location for placement of
 *     pasted annotation.
 */
void
UserInputModeAnnotations::pasteAnnotationFromAnnotationClipboard(const MouseEvent& mouseEvent)
{
    Annotation* newPastedAnnotation = AnnotationPasteDialog::pasteAnnotationOnClipboard(mouseEvent,
                                                                                        m_browserWindowIndex);
    if (newPastedAnnotation != NULL) {
        selectAnnotation(newPastedAnnotation);
        
        DisplayPropertiesAnnotation* dpa = GuiManager::get()->getBrain()->getDisplayPropertiesAnnotation();
        dpa->updateForNewAnnotation(newPastedAnnotation);
    }
    
    setMode(MODE_SELECT);
    
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    EventManager::get()->sendEvent(EventUserInterfaceUpdate().getPointer());
}

/**
 * Paste the annotation from the annotation clipboard and allow user to change its space.
 *
 * @param mouseEvent
 *     Mouse event containing XY location for placement of
 *     pasted annotation.
 */
void
UserInputModeAnnotations::pasteAnnotationFromAnnotationClipboardAndChangeSpace(const MouseEvent& mouseEvent)
{
    Annotation* newPastedAnnotation = AnnotationPasteDialog::pasteAnnotationOnClipboardChangeSpace(mouseEvent);
    if (newPastedAnnotation != NULL) {
        selectAnnotation(newPastedAnnotation);
        
        DisplayPropertiesAnnotation* dpa = GuiManager::get()->getBrain()->getDisplayPropertiesAnnotation();
        dpa->updateForNewAnnotation(newPastedAnnotation);
    }
    
    setMode(MODE_SELECT);
    
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
    EventManager::get()->sendEvent(EventUserInterfaceUpdate().getPointer());
}


/* =============================================================================== */

/**
 * Constructor.
 *
 * @param annotationFile
 *     File for annotation.
 * @param annotationSpace
 *     Space for annotation being created.
 * @param annotationType
 *     Type of annotation being created.
 * @param mouseEvent
 *     Mouse event.
 */
UserInputModeAnnotations::NewMouseDragCreateAnnotation::NewMouseDragCreateAnnotation(AnnotationFile* annotationFile,
                                                                                     const AnnotationCoordinateSpaceEnum::Enum annotationSpace,
                                                                                     const AnnotationTypeEnum::Enum annotationType,
                                                                                     const MouseEvent& mouseEvent)
{
    BrainOpenGLViewportContent* vpContent = mouseEvent.getViewportContent();
    CaretAssert(vpContent);
    if (vpContent == NULL) {
        CaretLogSevere("Viewport content is invalid.");
        return;
    }
    int32_t windowViewport[4];
    vpContent->getWindowViewport(windowViewport);
    m_windowOriginX = windowViewport[0];
    m_windowOriginY = windowViewport[1];
    
    m_mousePressWindowX = mouseEvent.getX() - m_windowOriginX;
    m_mousePressWindowY = mouseEvent.getY() - m_windowOriginY;
    
    int viewport[4];
    vpContent->getWindowViewport(viewport);
    m_windowWidth  = viewport[2];
    m_windowHeight = viewport[3];
    
    
    m_annotationFile = annotationFile;
    m_annotation = Annotation::newAnnotationOfType(annotationType,
                                                   AnnotationAttributesDefaultTypeEnum::USER);
    m_annotation->setCoordinateSpace(annotationSpace);
    CaretAssert(m_annotation);

    AnnotationMultiCoordinateShape* multiCoordShape = m_annotation->castToMultiCoordinateShape();
    AnnotationOneCoordinateShape* oneCoordShape     = m_annotation->castToOneCoordinateShape();
    AnnotationTwoCoordinateShape* twoCoordShape     = m_annotation->castToTwoCoordinateShape();

    Vector3D mouseCoord3D(mouseEvent.getPressedX(),
                          mouseEvent.getPressedY(),
                          0.0);
    if (twoCoordShape != NULL) {
        setCoordinate(twoCoordShape->getStartCoordinate(),
                      m_mousePressWindowX,
                      m_mousePressWindowY);
        m_drawingCoordinates.push_back(mouseCoord3D);

        setCoordinate(twoCoordShape->getEndCoordinate(),
                      m_mousePressWindowX,
                      m_mousePressWindowY);
        m_drawingCoordinates.push_back(mouseCoord3D);

        CaretAssert(m_drawingCoordinates.size() == 2);
    }
    else if (oneCoordShape != NULL) {
        setCoordinate(oneCoordShape->getCoordinate(),
                      m_mousePressWindowX,
                      m_mousePressWindowY);
        oneCoordShape->setWidth(1.0);
        oneCoordShape->setHeight(1.0);
        
        m_drawingCoordinates.push_back(mouseCoord3D);
        CaretAssert(m_drawingCoordinates.size() == 1);
    }
    else if (multiCoordShape != NULL) {
        AnnotationCoordinate* ac = new AnnotationCoordinate(AnnotationAttributesDefaultTypeEnum::USER);
        setCoordinate(ac, m_mousePressWindowX, m_mousePressWindowY);
        multiCoordShape->addCoordinate(ac);
        
        m_drawingCoordinates.push_back(mouseCoord3D);
    }
    else {
        CaretAssert(0);
    }
    
    
    if ((m_annotation->getLineColor() == CaretColorEnum::NONE)
        && (m_annotation->getBackgroundColor() == CaretColorEnum::NONE)) {
        m_annotation->setLineColor(CaretColorEnum::RED);
    }
}

/**
 * Destructor.
 */
UserInputModeAnnotations::NewMouseDragCreateAnnotation::~NewMouseDragCreateAnnotation()
{
    delete m_annotation;
}

/**
 * Update with current mouse location.
 *
 * @param mouseEvent
 *     Mouse event.
 * @param mouseWindowX
 *     Mouse window X-coordinate
 * @param mouseWindowY
 *     Mouse window Y-coordinate
 */
void
UserInputModeAnnotations::NewMouseDragCreateAnnotation::update(const MouseEvent& mouseEvent,
                                                               const int32_t mouseWindowXIn,
                                                               const int32_t mouseWindowYIn)
{
    int32_t mouseWindowX = mouseWindowXIn - m_windowOriginX;
    int32_t mouseWindowY = mouseWindowYIn - m_windowOriginY;
    
    AnnotationTwoCoordinateShape* twoCoordShape = m_annotation->castToTwoCoordinateShape();
    AnnotationOneCoordinateShape* oneCoordShape = m_annotation->castToOneCoordinateShape();
    AnnotationMultiCoordinateShape* multiCoordShape = m_annotation->castToMultiCoordinateShape();
    
    Vector3D mouseCoord3D(mouseEvent.getX(),
                          mouseEvent.getY(),
                          0.0);

    if (twoCoordShape != NULL) {
        setCoordinate(twoCoordShape->getEndCoordinate(),
                      mouseWindowX,
                      mouseWindowY);
        CaretAssertVectorIndex(m_drawingCoordinates, 1);
        m_drawingCoordinates[1] = mouseCoord3D;
    }
    else if (oneCoordShape != NULL) {
        const float minX = std::min(m_mousePressWindowX,
                                    mouseWindowX);
        const float maxX = std::max(m_mousePressWindowX,
                                    mouseWindowX);

        const float minY = std::min(m_mousePressWindowY,
                                    mouseWindowY);
        const float maxY = std::max(m_mousePressWindowY,
                                    mouseWindowY);

        const float x = (minX + maxX) / 2.0;
        const float y = (minY + maxY) / 2.0;
        
        const float width  = maxX - minX;
        const float height = maxY - minY;
        
        const float relativeWidth = 100.0 * ((m_windowWidth > 0.0)
                                             ? (width / static_cast<float>(m_windowWidth))
                                             : 0.01);
        const float relativeHeight = 100.0 * ((m_windowHeight > 0.0)
                                              ? (height / static_cast<float>(m_windowHeight))
                                              : 0.01);
        
        AnnotationCoordinate* coord = oneCoordShape->getCoordinate();
        setCoordinate(coord, x, y);
        oneCoordShape->setWidth(relativeWidth);
        oneCoordShape->setHeight(relativeHeight);
        
        CaretAssertVectorIndex(m_drawingCoordinates, 0);
        m_drawingCoordinates[0][0] = ((mouseEvent.getPressedX() + mouseEvent.getX()) / 2.0);
        m_drawingCoordinates[0][1] = ((mouseEvent.getPressedY() + mouseEvent.getY()) / 2.0);
    }
    else if (multiCoordShape != NULL) {
        AnnotationCoordinate* ac = new AnnotationCoordinate(AnnotationAttributesDefaultTypeEnum::USER);
        setCoordinate(ac, mouseWindowX, mouseWindowY);
        multiCoordShape->addCoordinate(ac);
        
        m_drawingCoordinates.push_back(mouseCoord3D);
    }
    else {
        CaretAssert(0);
    }
}

/**
 * Set the given coordinate to the XY-location.
 *
 * @param coordinate
 *     Coordinate that is set.
 * @param x
 *     Mouse window X-coordinate
 * @param y
 *     Mouse window Y-coordinate
 */
void
UserInputModeAnnotations::NewMouseDragCreateAnnotation::setCoordinate(AnnotationCoordinate* coordinate,
                                                                      const int32_t x,
                                                                      const int32_t y)
{
    const float relativeX = 100.0 * ((m_windowWidth > 0.0)
                                     ? (x / m_windowWidth)
                                     : 0.01);
    const float relativeY = 100.0 * ((m_windowHeight > 0.0)
                                     ? (y / m_windowHeight)
                                     : 0.01);
    
    coordinate->setXYZ(relativeX,
                       relativeY,
                       0.0);
}

/**
 * @return New annotation being drawn by the user.
 */
const Annotation*
UserInputModeAnnotations::NewMouseDragCreateAnnotation::getAnnotation() const
{
    return m_annotation;
}

/**
 * @return The drawing coordinates
 */
const std::vector<Vector3D>&
UserInputModeAnnotations::NewMouseDragCreateAnnotation::getDrawingCoordinates() const
{
    return m_drawingCoordinates;
}


