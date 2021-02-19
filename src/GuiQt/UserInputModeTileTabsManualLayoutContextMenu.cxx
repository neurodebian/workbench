
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

#define __USER_INPUT_MODE_TILE_TABS_MANUAL_LAYOUT_CONTEXT_MENU_DECLARE__
#include "UserInputModeTileTabsManualLayoutContextMenu.h"
#undef __USER_INPUT_MODE_TILE_TABS_MANUAL_LAYOUT_CONTEXT_MENU_DECLARE__

#include <cmath>

#include <QLineEdit>

#include "AnnotationBrowserTab.h"
#include "AnnotationManager.h"
#include "AnnotationRedoUndoCommand.h"
#include "Brain.h"
#include "BrainOpenGLWidget.h"
#include "BrainBrowserWindow.h"
#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "EventBrowserTabGetAll.h"
#include "EventBrowserWindowTileTabOperation.h"
#include "EventGraphicsUpdateAllWindows.h"
#include "EventManager.h"
#include "EventUserInterfaceUpdate.h"
#include "GuiManager.h"
#include "MathFunctions.h"
#include "MouseEvent.h"
#include "UserInputModeTileTabsManualLayout.h"
#include "WuQDataEntryDialog.h"
#include "WuQMessageBox.h"

using namespace caret;


    
/**
 * \class caret::UserInputModeTileTabsManualLayoutContextMenu
 * \brief Context (pop-up) menu for User Input Manual Tile Tabs Layout.
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param userInputTileTabsManualLayout
 *    The manual layout input processor
 * @param mouseEvent
 *    The mouse event that caused display of this menu.
 * @param browserTabContent
 *    Content of browser tab.
 * @param parentOpenGLWidget
 *    Parent OpenGL Widget on which the menu is displayed.
 */
UserInputModeTileTabsManualLayoutContextMenu::UserInputModeTileTabsManualLayoutContextMenu(UserInputModeTileTabsManualLayout* userInputTileTabsManualLayout,
                                                                         const MouseEvent& mouseEvent,
                                                                         BrowserTabContent* browserTabContent,
                                                                         BrainOpenGLWidget* parentOpenGLWidget)
: QMenu(parentOpenGLWidget),
m_userInputTileTabsManualLayout(userInputTileTabsManualLayout),
m_mouseEvent(mouseEvent),
m_browserTabContent(browserTabContent),
m_parentOpenGLWidget(parentOpenGLWidget)
{
    CaretAssert(m_userInputTileTabsManualLayout);
    
    const int32_t browserWindexIndex = m_mouseEvent.getBrowserWindowIndex();
    AnnotationManager* annotationManager = GuiManager::get()->getBrain()->getAnnotationManager();
    std::vector<Annotation*> selectedAnnotations = annotationManager->getAnnotationsSelectedForEditing(browserWindexIndex);
    for (auto ann : selectedAnnotations) {
        AnnotationBrowserTab* abt = dynamic_cast<AnnotationBrowserTab*>(ann);
        if (abt != NULL) {
            m_selectedBrowserTabAnnotations.push_back(abt);
        }
    }
    const bool oneTabSelectedFlag((m_selectedBrowserTabAnnotations.size() == 1)
                                  && (m_browserTabContent != NULL));
    
    /*
     * De/Select All tabs
     */
    QAction* deselectAction = addAction(BrainBrowserWindowEditMenuItemEnum::toGuiName(BrainBrowserWindowEditMenuItemEnum::DESELECT_ALL),
                                        this, SLOT(deselectAllAnnotations()));
    deselectAction->setEnabled( ! m_selectedBrowserTabAnnotations.empty());
    addAction(BrainBrowserWindowEditMenuItemEnum::toGuiName(BrainBrowserWindowEditMenuItemEnum::SELECT_ALL),
              this, SLOT(selectAllAnnotations()));
    
    
    /*
     * Separator
     */
    addSeparator();
    
    /*
     * Order Operations
     */
    QAction* bringToFrontAction = addAction("Bring to Front");
    QObject::connect(bringToFrontAction, &QAction::triggered,
                     this, &UserInputModeTileTabsManualLayoutContextMenu::applyOrderBringToFront);
    QAction* bringForwardAction = addAction("Bring Forward");
    QObject::connect(bringForwardAction, &QAction::triggered,
                     this, &UserInputModeTileTabsManualLayoutContextMenu::applyOrderBringForward);
    QAction* sendToBackAction = addAction("Send to Back");
    QObject::connect(sendToBackAction, &QAction::triggered,
                     this, &UserInputModeTileTabsManualLayoutContextMenu::applyOrderSendToBack);
    QAction* sendBackwardAction = addAction("Send Backward");
    QObject::connect(sendBackwardAction, &QAction::triggered,
                     this, &UserInputModeTileTabsManualLayoutContextMenu::applyOrderSendBackward);

    bringToFrontAction->setEnabled(oneTabSelectedFlag);
    bringForwardAction->setEnabled(oneTabSelectedFlag);
    sendToBackAction->setEnabled(oneTabSelectedFlag);
    sendBackwardAction->setEnabled(oneTabSelectedFlag);
    
    const bool showGroupingOptionsFlag(false);
    if (showGroupingOptionsFlag) {
        addSeparator();
        
        /*
         * Group annotations
         */
        QAction* groupAction = addAction(AnnotationGroupingModeEnum::toGuiName(AnnotationGroupingModeEnum::GROUP),
                                         this, SLOT(applyGroupingGroup()));
        groupAction->setEnabled(annotationManager->isGroupingModeValid(browserWindexIndex,
                                                                       AnnotationGroupingModeEnum::GROUP));
        
        /*
         * Ungroup annotations
         */
        QAction* ungroupAction = addAction(AnnotationGroupingModeEnum::toGuiName(AnnotationGroupingModeEnum::UNGROUP),
                                           this, SLOT(applyGroupingUngroup()));
        ungroupAction->setEnabled(annotationManager->isGroupingModeValid(browserWindexIndex,
                                                                         AnnotationGroupingModeEnum::UNGROUP));
        
        /*
         * Regroup annotations
         */
        QAction* regroupAction = addAction(AnnotationGroupingModeEnum::toGuiName(AnnotationGroupingModeEnum::REGROUP),
                                           this, SLOT(applyGroupingRegroup()));
        regroupAction->setEnabled(annotationManager->isGroupingModeValid(browserWindexIndex,
                                                                         AnnotationGroupingModeEnum::REGROUP));
    }

    addSeparator();
    
    /*
     * Insert new tab option
     */
    QAction* insertNewTabAction = addAction("Insert New Tab");
    QObject::connect(insertNewTabAction, &QAction::triggered,
                     this, &UserInputModeTileTabsManualLayoutContextMenu::processInsertNewTabMenuItem);
    
    addSeparator();
    
    /*
     * Expand option
     */
    QAction* expendToFillAction = addAction(getShinkAndExpandTabMenuItemText());
    QObject::connect(expendToFillAction, &QAction::triggered,
                     this, &UserInputModeTileTabsManualLayoutContextMenu::processShrinkAndExpandTabMenuItem);
    expendToFillAction->setEnabled(oneTabSelectedFlag);
}

/**
 * Destructor.
 */
UserInputModeTileTabsManualLayoutContextMenu::~UserInputModeTileTabsManualLayoutContextMenu()
{
}

/**
 * @return Text for "shrink and expand tab" menu item (this function is available from
 * this context menu and the Arrange Menu in the toolbar.  This function ensures the
 * same text in both menus.
 */
AString
UserInputModeTileTabsManualLayoutContextMenu::getShinkAndExpandTabMenuItemText()
{
    return "Move/Resize Tab to Fill Empty Space";
}

/**
 * Deselect all annotations in the window.
 */
void
UserInputModeTileTabsManualLayoutContextMenu::deselectAllAnnotations()
{
    m_userInputTileTabsManualLayout->processDeselectAllAnnotations();
}

/**
 * Select all annotations in the window.
 */
void
UserInputModeTileTabsManualLayoutContextMenu::selectAllAnnotations()
{
    m_userInputTileTabsManualLayout->processSelectAllAnnotations();
}

/**
 * Group annotations.
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyGroupingGroup()
{
    applyGrouping(AnnotationGroupingModeEnum::GROUP);
}

/**
 * Ungroup annotations.
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyGroupingRegroup()
{
    applyGrouping(AnnotationGroupingModeEnum::REGROUP);
}

/**
 * Regroup annotations.
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyGroupingUngroup()
{
    applyGrouping(AnnotationGroupingModeEnum::UNGROUP);
}

/**
 * Apply grouping selection.
 *
 * @param grouping
 *     Selected grouping.
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyGrouping(const AnnotationGroupingModeEnum::Enum grouping)
{
    AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
    
    AString errorMessage;
    if ( ! annMan->applyGroupingMode(m_userInputTileTabsManualLayout->getUserInputMode(),
                                     m_mouseEvent.getBrowserWindowIndex(),
                                     grouping,
                                     errorMessage)) {
        WuQMessageBox::errorOk(this,
                               errorMessage);
    }
    
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());
}

/**
 * Called when expand to fill space menu item is selected
 */
void
UserInputModeTileTabsManualLayoutContextMenu::processShrinkAndExpandTabMenuItem()
{
    const int32_t browserWindowIndex = m_mouseEvent.getBrowserWindowIndex();
    BrainBrowserWindow* window = GuiManager::get()->getBrowserWindowByWindowIndex(browserWindowIndex);
    CaretAssert(window);
    std::vector<BrowserTabContent*> allTabContent;
    window->getAllTabContent(allTabContent);
    
    AnnotationManager* annMan = GuiManager::get()->getBrain()->getAnnotationManager();
    AString errorMessage;
    const bool result = annMan->shrinkAndExpandSelectedBrowserTabAnnotation(allTabContent,
                                                                   browserWindowIndex,
                                                                   m_userInputTileTabsManualLayout->getUserInputMode(),
                                                                   errorMessage);
    EventManager::get()->sendSimpleEvent(EventTypeEnum::EVENT_ANNOTATION_TOOLBAR_UPDATE);
    EventManager::get()->sendEvent(EventGraphicsUpdateAllWindows().getPointer());

    if ( ! result) {
        WuQMessageBox::errorOk(this,
                               errorMessage);
    }    
}

/**
 * Called when bring to front is selected
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyOrderBringToFront()
{
    processWindowTileTabOperation(EventBrowserWindowTileTabOperation::OPERATION_ORDER_BRING_TO_FRONT);
}

/**
 * Called when bring forward is selected
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyOrderBringForward()
{
    processWindowTileTabOperation(EventBrowserWindowTileTabOperation::OPERATION_ORDER_BRING_FORWARD);
}

/**
 * Called when send to back is selected
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyOrderSendToBack()
{
    processWindowTileTabOperation(EventBrowserWindowTileTabOperation::OPERATION_ORDER_SEND_TO_BACK);
}

/**
 * Called when send backward is selected
 */
void
UserInputModeTileTabsManualLayoutContextMenu::applyOrderSendBackward()
{
    processWindowTileTabOperation(EventBrowserWindowTileTabOperation::OPERATION_ORDER_SEND_BACKWARD);
}

/**
 * Insert a new tab at mouse location
 */
void
UserInputModeTileTabsManualLayoutContextMenu::processInsertNewTabMenuItem()
{
    processWindowTileTabOperation(EventBrowserWindowTileTabOperation::OPERATION_MANUAL_NEW_TAB);
}


/**
 * Called to process a tile tab operation
 *
 * @param operation
 *     The operation
 */
void
UserInputModeTileTabsManualLayoutContextMenu::processWindowTileTabOperation(const EventBrowserWindowTileTabOperation::Operation operation)
{
    switch (operation) {
        case EventBrowserWindowTileTabOperation::OPERATION_GRID_NEW_TAB_AFTER:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_GRID_NEW_TAB_BEFORE:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_MANUAL_NEW_TAB:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_ORDER_BRING_FORWARD:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_ORDER_BRING_TO_FRONT:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_ORDER_SEND_BACKWARD:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_ORDER_SEND_TO_BACK:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_REPLACE_TABS:
            break;
        case EventBrowserWindowTileTabOperation::OPERATION_SELECT_TAB:
            break;
    }
    
    std::vector<BrowserTabContent*> emptyBrowserTabs;
    
    int tabIndex(-1);
    if (m_browserTabContent != NULL) {
        tabIndex = m_browserTabContent->getTabNumber();
    }
    int windowViewport[4] { 0, 0, m_parentOpenGLWidget->width(), m_parentOpenGLWidget->height() };
    EventBrowserWindowTileTabOperation tileTabOperation(operation,
                                                        m_parentOpenGLWidget,
                                                        m_mouseEvent.getBrowserWindowIndex(),
                                                        tabIndex,
                                                        windowViewport,
                                                        m_mouseEvent.getX(),
                                                        m_mouseEvent.getY(),
                                                        emptyBrowserTabs);
    EventManager::get()->sendEvent(tileTabOperation.getPointer());
}



