
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

#define __EVENT_BROWSER_WINDOW_TILE_TAB_OPERATION_DECLARE__
#include "EventBrowserWindowTileTabOperation.h"
#undef __EVENT_BROWSER_WINDOW_TILE_TAB_OPERATION_DECLARE__

#include <QWidget>

#include "CaretAssert.h"
#include "EventManager.h"
#include "EventTypeEnum.h"
#include "WuQMessageBox.h"
using namespace caret;


    
/**
 * \class caret::EventBrowserWindowTileTabOperation 
 * \brief Operations for tab changes
 * \ingroup GuiQt
 */

/**
 * Constructor.
 *
 * @param operation
 *    The operation
 * @param parentWidget
 *    Parent widget for error dialogs.
 * @param windowIndex
 *    Index of the window.
 * @param browserTabIndex
 *    Index of the browser tab.
 * @param windowViewport
 *    The window viewport
 * @param mouseX
 *    X-position of mouse
 * @param mouseY
 *    Y-position of mouse
 * @param browserTabsForReplaceOperation
 *    Tabs for replacement
 */
EventBrowserWindowTileTabOperation::EventBrowserWindowTileTabOperation(const Operation operation,
                                                                       QWidget* parentWidget,
                                                                       const int32_t windowIndex,
                                                                       const int32_t browserTabIndex,
                                                                       const int32_t windowViewport[4],
                                                                       const int32_t mouseX,
                                                                       const int32_t mouseY,
                                                                       const std::vector<BrowserTabContent*>& browserTabsForReplaceOperation)
: Event(EventTypeEnum::EVENT_BROWSER_WINDOW_TILE_TAB_OPERATION),
m_operation(operation),
m_parentWidget(parentWidget),
m_windowIndex(windowIndex),
m_browserTabIndex(browserTabIndex),
m_mouseX(mouseX),
m_mouseY(mouseY),
m_browserTabsForReplaceOperation(browserTabsForReplaceOperation)
{
    CaretAssert(m_parentWidget);
    CaretAssert(m_windowIndex >= 0);

    m_windowViewport[0] = windowViewport[0];
    m_windowViewport[1] = windowViewport[1];
    m_windowViewport[2] = windowViewport[2];
    m_windowViewport[3] = windowViewport[3];

    switch (m_operation) {
        case OPERATION_GRID_NEW_TAB_AFTER:
            CaretAssert(m_browserTabIndex >= 0);
            break;
        case OPERATION_GRID_NEW_TAB_BEFORE:
            CaretAssert(m_browserTabIndex >= 0);
            break;
        case OPERATION_MANUAL_NEW_TAB:
            break;
        case OPERATION_ORDER_BRING_TO_FRONT:
            CaretAssert(m_browserTabIndex >= 0);
            break;
        case OPERATION_ORDER_BRING_FORWARD:
            CaretAssert(m_browserTabIndex >= 0);
            break;
        case OPERATION_ORDER_SEND_TO_BACK:
            CaretAssert(m_browserTabIndex >= 0);
            break;
        case OPERATION_ORDER_SEND_BACKWARD:
            CaretAssert(m_browserTabIndex >= 0);
            break;
        case OPERATION_REPLACE_TABS:
            break;
        case OPERATION_SELECT_TAB:
            CaretAssert(m_browserTabIndex >= 0);
            break;
    }
}

/**
 * Destructor.
 */
EventBrowserWindowTileTabOperation::~EventBrowserWindowTileTabOperation()
{
}

/**
 * Operation that selects the given tab in the given window.
 *
 * @param parentWidget
 *    Parent widget for error dialogs.
 * @param windowIndex
 *    Index of the window.
 * @param browserTabIndex
 *    Index of the browser tab.
 */
void
EventBrowserWindowTileTabOperation::selectTabInWindow(QWidget* parentWidget,
                                                      const int32_t windowIndex,
                                                      const int32_t browserTabIndex)
{
    const int32_t dummyMouseX(-1);
    const int32_t dummyMouseY(-1);
    const int32_t dummyWindowViewport[4] { -1, -1, -1, -1 };
    std::vector<BrowserTabContent*> emptyBrowserTabs;
    EventBrowserWindowTileTabOperation tabOperation(Operation::OPERATION_SELECT_TAB,
                                                    parentWidget,
                                                    windowIndex,
                                                    browserTabIndex,
                                                    dummyWindowViewport,
                                                    dummyMouseX,
                                                    dummyMouseY,
                                                    emptyBrowserTabs);
    
    EventManager::get()->sendEvent(tabOperation.getPointer());
    
    if (tabOperation.getEventProcessCount() <= 0) {
        WuQMessageBox::errorOk(parentWidget,
                               "Tab not created, invalid window or tab index");
    }
}

/**
 * @return The mode
 */
EventBrowserWindowTileTabOperation::Operation
EventBrowserWindowTileTabOperation::getOperation() const
{
    return m_operation;
}

/**
 * @return The window index.
 */
int32_t
EventBrowserWindowTileTabOperation::getWindowIndex() const
{
    return m_windowIndex;
}

/**
 * @return The browser tab index.
 */
int32_t
EventBrowserWindowTileTabOperation::getBrowserTabIndex() const
{
    return m_browserTabIndex;
}

/**
 * @return Get the browser tabs for a replace tabs operation.
 */
const std::vector<BrowserTabContent*>
EventBrowserWindowTileTabOperation::getBrowserTabsForReplaceOperation() const
{
    return m_browserTabsForReplaceOperation;
}

/**
 * @return Mouse X-coordinate (invalid if negative)
 */
int
EventBrowserWindowTileTabOperation::getMouseX() const
{
    return m_mouseX;
}

/**
 * @return Mouse Y-coordinate (invalid if negative)
 */
int
EventBrowserWindowTileTabOperation::getMouseY() const
{
    return m_mouseY;
}

/**
 * Get the window viewport
 *
 * @param windowViewport
 *     Output containing window viewport (negative values if invalid)
 */
void
EventBrowserWindowTileTabOperation::getWindowViewport(int32_t windowViewportOut[4]) const
{
    windowViewportOut[0] = m_windowViewport[0];
    windowViewportOut[1] = m_windowViewport[1];
    windowViewportOut[2] = m_windowViewport[2];
    windowViewportOut[3] = m_windowViewport[3];
}

