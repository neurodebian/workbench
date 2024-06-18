
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

#define __BRAIN_OPEN_G_L_VIEWPORT_CONTENT_DECLARE__
#include "BrainOpenGLViewportContent.h"
#undef __BRAIN_OPEN_G_L_VIEWPORT_CONTENT_DECLARE__

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "AnnotationBrowserTab.h"
#include "BrowserTabContent.h"
#include "BrowserWindowContent.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "EventBrowserWindowContent.h"
#include "EventManager.h"
#include "EventSpacerTabGet.h"
#include "GapsAndMargins.h"
#include "GraphicsObjectToWindowTransform.h"
#include "MathFunctions.h"
#include "MediaFile.h"
#include "MediaOverlay.h"
#include "MediaOverlaySet.h"
#include "ModelSurfaceMontage.h"
#include "SpacerTabContent.h"
#include "SurfaceMontageConfigurationAbstract.h"
#include "TileTabsBrowserTabGeometry.h"
#include "TileTabsLayoutGridConfiguration.h"
#include "TileTabsLayoutManualConfiguration.h"

using namespace caret;


    
/**
 * \class BrainOpenGLViewportContent 
 * \brief Dimensions and model for a viewport in the graphics window.
 *
 * Dimensions and model for a viewport in the graphics window.
 */


/**
 * Constructor.
 *
 * @param windowBeforeAspectLockingViewport
 *    Viewport of window before aspect locking is applied
 * @param windowViewport
 *    Viewport of WINDOW in which drawing takes place.
 * @param tabViewportManualLayoutBeforeAspectLocking
 *    Viewport before locking in a manual layout (for grid, same as tab viewport)
 * @param tabViewport
 *    Viewport for TAB in which drawing takes place.
 * @param modelViewport
 *    Viewport for MODEL in which drawing takes place.
 * @param windowIndex
 *    Index of browser window.
 * @param highlightTabFlag
 *    True indicates that the tab is highlighted (used in
 *    Tile Tabs mode so user knows graphics region corresponding
 *    to the tab that was recently selected).
 * @param browserTabContent
 *    Browser Tab content that is being drawn (if not NULL)
 * @param spacerTabContent
 *    Spacer Tab content that is being drawn (if not NULL)
 */
BrainOpenGLViewportContent::BrainOpenGLViewportContent(const int32_t windowBeforeAspectLockingViewport[4],
                                                       const int windowViewport[4],
                                                       const int tabViewportManualLayoutBeforeAspectLocking[4],
                                                       const int tabViewport[4],
                                                       const int modelViewport[4],
                                                       const int windowIndex,
                                                       const bool highlightTabFlag,
                                                       BrowserTabContent* browserTabContent,
                                                       SpacerTabContent* spacerTabContent)
: CaretObject(),
m_windowIndex(windowIndex),
m_highlightTab(highlightTabFlag),
m_browserTabContent(browserTabContent),
m_spacerTabContent(spacerTabContent),
m_windowBeforeAspectLockingX(windowBeforeAspectLockingViewport[0]),
m_windowBeforeAspectLockingY(windowBeforeAspectLockingViewport[1]),
m_windowBeforeAspectLockingWidth(windowBeforeAspectLockingViewport[2]),
m_windowBeforeAspectLockingHeight(windowBeforeAspectLockingViewport[3])
{
    /*
     * Note all members initialized in the header file
     */
    m_windowX      = windowViewport[0];
    m_windowY      = windowViewport[1];
    m_windowWidth  = windowViewport[2];
    m_windowHeight = windowViewport[3];
    
    m_tabX      = tabViewport[0];
    m_tabY      = tabViewport[1];
    m_tabWidth  = tabViewport[2];
    m_tabHeight = tabViewport[3];
    
    for (int32_t i = 0; i < 4; i++) {
        m_tabViewportManualLayoutBeforeAspectLocking[i] = tabViewportManualLayoutBeforeAspectLocking[i];
    }
    
    m_modelX      = modelViewport[0];
    m_modelY      = modelViewport[1];
    m_modelWidth  = modelViewport[2];
    m_modelHeight = modelViewport[3];
}

/**
 * Destructor.
 */
BrainOpenGLViewportContent::~BrainOpenGLViewportContent()
{
    
}

/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
BrainOpenGLViewportContent::BrainOpenGLViewportContent(const BrainOpenGLViewportContent& obj)
: CaretObject(obj),
m_windowIndex(obj.m_windowIndex),
m_highlightTab(obj.m_highlightTab),
m_windowBeforeAspectLockingX(obj.m_windowBeforeAspectLockingX),
m_windowBeforeAspectLockingY(obj.m_windowBeforeAspectLockingY),
m_windowBeforeAspectLockingWidth(obj.m_windowBeforeAspectLockingWidth),
m_windowBeforeAspectLockingHeight(obj.m_windowBeforeAspectLockingHeight)
{
    /*
     * Note all members initialized in the header file
     */
    this->copyHelperBrainOpenGLViewportContent(obj);
}

/**
 * Assignment operator.
 * @param obj
 *    Data copied from obj to this.
 * @return
 *    Reference to this object.
 */
BrainOpenGLViewportContent&
BrainOpenGLViewportContent::operator=(const BrainOpenGLViewportContent& obj)
{
    if (this != &obj) {
        CaretObject::operator=(obj);
        this->copyHelperBrainOpenGLViewportContent(obj);
    }
    return *this;
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void
BrainOpenGLViewportContent::copyHelperBrainOpenGLViewportContent(const BrainOpenGLViewportContent& obj)
{
    m_chartDataProjectionMatrix = obj.m_chartDataProjectionMatrix;
    m_chartDataModelViewMatrix  = obj.m_chartDataModelViewMatrix;
    m_chartDataX      = obj.m_chartDataX;
    m_chartDataY      = obj.m_chartDataY;
    m_chartDataWidth  = obj.m_chartDataWidth;
    m_chartDataHeight = obj.m_chartDataHeight;
    m_chartDataViewportValidFlag = obj.m_chartDataViewportValidFlag;
    m_modelX       = obj.m_modelX;
    m_modelY       = obj.m_modelY;
    m_modelWidth   = obj.m_modelWidth;
    m_modelHeight  = obj.m_modelHeight;
    m_tabX         = obj.m_tabX;
    m_tabY         = obj.m_tabY;
    m_tabWidth     = obj.m_tabWidth;
    m_tabHeight    = obj.m_tabHeight;
    m_windowX      = obj.m_windowX;
    m_windowY      = obj.m_windowY;
    m_windowWidth  = obj.m_windowWidth;
    m_windowHeight = obj.m_windowHeight;
    m_windowBeforeAspectLockingX      = obj.m_windowBeforeAspectLockingX;
    m_windowBeforeAspectLockingY      = obj.m_windowBeforeAspectLockingY;
    m_windowBeforeAspectLockingWidth  = obj.m_windowBeforeAspectLockingWidth;
    m_windowBeforeAspectLockingHeight = obj.m_windowBeforeAspectLockingHeight;

    for (int32_t i = 0; i < 4; i++) {
        m_tabViewportManualLayoutBeforeAspectLocking[i] = obj.m_tabViewportManualLayoutBeforeAspectLocking[i];
    }

    m_browserTabContent = obj.m_browserTabContent;
    m_spacerTabContent  = obj.m_spacerTabContent;
    m_histologyGraphicsObjectToWindowTransform.reset();
    if (obj.m_histologyGraphicsObjectToWindowTransform) {
        m_histologyGraphicsObjectToWindowTransform.reset(new GraphicsObjectToWindowTransform(*obj.m_histologyGraphicsObjectToWindowTransform));
    }
    m_mediaGraphicsObjectToWindowTransform.reset();
    if (obj.m_mediaGraphicsObjectToWindowTransform) {
        m_mediaGraphicsObjectToWindowTransform.reset(new GraphicsObjectToWindowTransform(*obj.m_mediaGraphicsObjectToWindowTransform));
    }
    m_volumeAxialGraphicsObjectToWindowTransform.reset();
    if (obj.m_volumeAxialGraphicsObjectToWindowTransform) {
        m_volumeAxialGraphicsObjectToWindowTransform.reset(new GraphicsObjectToWindowTransform(*obj.m_volumeAxialGraphicsObjectToWindowTransform));
    }
    m_volumeCoronalGraphicsObjectToWindowTransform.reset();
    if (obj.m_volumeCoronalGraphicsObjectToWindowTransform) {
        m_volumeCoronalGraphicsObjectToWindowTransform.reset(new GraphicsObjectToWindowTransform(*obj.m_volumeCoronalGraphicsObjectToWindowTransform));
    }
    m_volumeParasagittalGraphicsObjectToWindowTransform.reset();
    if (obj.m_volumeParasagittalGraphicsObjectToWindowTransform) {
        m_volumeParasagittalGraphicsObjectToWindowTransform.reset(new GraphicsObjectToWindowTransform(*obj.m_volumeParasagittalGraphicsObjectToWindowTransform));
    }
}

/**
 * Adjust the width/height using the aspect ratio
 */
void
BrainOpenGLViewportContent::adjustWidthHeightForAspectRatio(const float aspectRatio,
                                                            int32_t& width,
                                                            int32_t& height)
{
    if (aspectRatio > 0.0) {
        int32_t heightInt = height;
        float widthFloat  = width;
        float heightFloat = height;
        
        float preferredHeightFloat = MathFunctions::round(widthFloat * aspectRatio);
        const int32_t preferredHeightInt = static_cast<int32_t>(preferredHeightFloat);
        if (heightInt == preferredHeightInt) {
            /*
             * Viewport matches aspect ratio so do not need to 
             * adjust the width and height.
             *
             * Due to floating point error, when lock is enabled,
             * the preferred height may be a very small difference
             * from the current height.  So rounding and then
             * converting to an int prevents the graphics region
             * from a small resizing.
             */
        }
        else if (preferredHeightFloat > heightFloat) {
            /*
             * Shrink width
             */
            const float percentage = heightFloat / preferredHeightFloat;
            width = static_cast<int32_t>(widthFloat * percentage);
        }
        else {
            /*
             * Shrink height
             */
            height = preferredHeightInt;
        }
    }
}


/**
 * Adjust the given viewport by applying the given aspect ratio.
 *
 * Sets the new height to be width * aspect ratio.  If this new height
 * is too tall, the viewport width and height is scaled down so that 
 * the height fits the original viewport size and the viewport is
 * horizontally centered.  If the new height is less than the original
 * height, the viewport is centered vertically.
 *
 * @param viewport
 *     The viewport
 * @param aspectRatio
 *     The aspect ratio (height ~= width * aspect ratio)
 */
void
BrainOpenGLViewportContent::adjustViewportForAspectRatio(int viewport[4],
                                                         const float aspectRatio)
{
    int32_t heightInt = viewport[3];
    float widthFloat  = viewport[2];
    float heightFloat = viewport[3];
    
    float preferredHeightFloat = MathFunctions::round(widthFloat * aspectRatio);
    const int32_t preferredHeightInt = static_cast<int32_t>(preferredHeightFloat);
    if (heightInt == preferredHeightInt) {
        /*
         * Due to floating point error, when lock is enabled,
         * the preferred height may be a very small difference
         * from the current height.  So rounding and then
         * converting to an int prevents the graphics region
         * from a small resizing.
         */
    }
    else if (preferredHeightFloat > heightFloat) {
        const float percentage = heightFloat / preferredHeightFloat;
        widthFloat *= percentage;
        
        const float xOffset = (viewport[2] - widthFloat) / 2.0;
        viewport[0] += xOffset;
        viewport[2]  = widthFloat;
    }
    else {
        const float yOffset = (viewport[3] - preferredHeightFloat) / 2.0;
        viewport[1] += yOffset;
        viewport[3] = preferredHeightFloat;
    }
}


/**
 * @return  True indicates that the tab is highlighted (used in
 *    Tile Tabs mode so user knows graphics region corresponding
 *    to the tab that was recently selected).
 */
bool
BrainOpenGLViewportContent::isTabHighlighted() const
{
    return m_highlightTab;
}

/**
 * Get the data bounds and viewport for drawing the chart data.
 *
 * @param chartDataProjectionMatrixOut
 *    Output into which projection matrix for drawing chart data loaded.
 * @param chartDataModelViewMatrixOut
 *    Output into which model viewing matrix for drawing chart data loaded.
 * @param chartViewportOut
 *    Output into which viewport is loaded.
 *    (x, y, width, height)
 * @return 
 *    True if the chart data viewport is valid.
 */
bool
BrainOpenGLViewportContent::getChartDataMatricesAndViewport(Matrix4x4& chartDataProjectionMatrixOut,
                                                         Matrix4x4& chartDataModelViewMatrixOut,
                                                         int chartViewportOut[4]) const
{
    chartDataProjectionMatrixOut = m_chartDataProjectionMatrix;
    chartDataModelViewMatrixOut  = m_chartDataModelViewMatrix;
    chartViewportOut[0] = m_chartDataX;
    chartViewportOut[1] = m_chartDataY;
    chartViewportOut[2] = m_chartDataWidth;
    chartViewportOut[3] = m_chartDataHeight;
    return m_chartDataViewportValidFlag;
}

/**
 * Set the viewport for drawing the chart data.
 *
 * @param chartDataProjectionMatrix
 *    Pojection matrix for drawing chart data.
 * @param chartDataModelViewMatrix
 *    Viewing matrix for drawing chart data.
 * @param chartDataViewport
 *    Viewport (x, y, width, height).
 */
void
BrainOpenGLViewportContent::setChartDataMatricesAndViewport(const Matrix4x4& chartDataProjectionMatrix,
                                                         const Matrix4x4& chartDataModelViewMatrix,
                                                         const int chartViewport[4]) const
{
    m_chartDataProjectionMatrix = chartDataProjectionMatrix;
    m_chartDataModelViewMatrix  = chartDataModelViewMatrix;
    m_chartDataX = chartViewport[0];
    m_chartDataY = chartViewport[1];
    m_chartDataWidth = chartViewport[2];
    m_chartDataHeight = chartViewport[3];
    m_chartDataViewportValidFlag = true;
}


/**
 * Get the viewport for drawing the model (has been reduced
 * from tab viewport by applying the margin).
 *
 * @param modelViewport
 *    Output into which model viewport dimensions are loaded.
 *    (x, y, width, height)
 */
void
BrainOpenGLViewportContent::getModelViewport(int modelViewport[4]) const
{
    modelViewport[0] = m_modelX;
    modelViewport[1] = m_modelY;
    modelViewport[2] = m_modelWidth;
    modelViewport[3] = m_modelHeight;
}

/**
 * Get the viewport for drawing the tab (includes margin).
 *
 * @param tabViewport
 *    Output into which tab viewport dimensions are loaded.
 *    (x, y, width, height)
 */
void
BrainOpenGLViewportContent::getTabViewportBeforeApplyingMargins(int tabViewportOut[4]) const
{
    tabViewportOut[0] = m_tabX;
    tabViewportOut[1] = m_tabY;
    tabViewportOut[2] = m_tabWidth;
    tabViewportOut[3] = m_tabHeight;
}
/**
 * Get the viewport for drawing a tab in a manual layout before application of aspect locking
 *
 * @param tabViewportOut
 *    Output into which tab viewport dimensions are loaded.
 *    (x, y, width, height)
 */
void
BrainOpenGLViewportContent::getTabViewportManualLayoutBeforeAspectLocking(int tabViewportOut[4]) const
{
    for (int32_t i = 0; i < 4; i++) {
        tabViewportOut[i] = m_tabViewportManualLayoutBeforeAspectLocking[i];
    }
}

/**
 * @param windowViewportOut
 *    Output into which window viewport dimensions are loaded.
 *    (x, y, width, height).  Viewport is AFTER aspect locking.
 */
void
BrainOpenGLViewportContent::getWindowViewport(int windowViewportOut[4]) const
{
    windowViewportOut[0] = m_windowX;
    windowViewportOut[1] = m_windowY;
    windowViewportOut[2] = m_windowWidth;
    windowViewportOut[3] = m_windowHeight;
}

/**
 * @param windowBeforeAspectLockingViewportOut
 *    Output into which window viewport dimensions are loaded.
 *    (x, y, width, height).  Viewport is BEFORE aspect locking.
 */
void
BrainOpenGLViewportContent::getWindowBeforeAspectLockingViewport(int windowBeforeAspectLockingViewportOut[4]) const
{
    windowBeforeAspectLockingViewportOut[0] = m_windowBeforeAspectLockingX;
    windowBeforeAspectLockingViewportOut[1] = m_windowBeforeAspectLockingY;
    windowBeforeAspectLockingViewportOut[2] = m_windowBeforeAspectLockingWidth;
    windowBeforeAspectLockingViewportOut[3] = m_windowBeforeAspectLockingHeight;
}

/**
 * @return The window index.
 */
int
BrainOpenGLViewportContent::getWindowIndex() const
{
    return m_windowIndex;
}

/**
 * @return Pointer to browser tab content in viewport (NULL if no browser tab in viewport)
 */
BrowserTabContent* 
BrainOpenGLViewportContent::getBrowserTabContent() const
{
    return m_browserTabContent;
}

/**
 * @return Pointer to spacer tab content in viewport (NULL if no spacer tab in viewport)
 */
SpacerTabContent*
BrainOpenGLViewportContent::getSpacerTabContent() const
{
    return m_spacerTabContent;
}

/**
 * @return Index of browser tab or -1 if there is not browser tab for this viewport.
 */
int32_t
BrainOpenGLViewportContent::getTabIndex() const
{
    if (m_browserTabContent != NULL) {
        return m_browserTabContent->getTabNumber();
    }
    
    return -1;
}


/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
BrainOpenGLViewportContent::toString() const
{
    const QString windowBeforeLockMsg = QString("   Window Before Lock x=%1 y=%2 w=%3 h=%4").arg(m_windowBeforeAspectLockingX).arg(m_windowBeforeAspectLockingY).arg(m_windowBeforeAspectLockingWidth).arg(m_windowBeforeAspectLockingHeight);
    const QString windowMsg = QString("   Window x=%1 y=%2 w=%3 h=%4").arg(m_windowX).arg(m_windowY).arg(m_windowWidth).arg(m_windowHeight);
    const QString tabMsg    = QString("   Tab    x=%1 y=%2 w=%3 h=%4").arg(m_tabX).arg(m_tabY).arg(m_tabWidth).arg(m_tabHeight);
    const QString modelMsg  = QString("   Model  x=%1 y=%2 w=%3 h=%4").arg(m_modelX).arg(m_modelY).arg(m_modelWidth).arg(m_modelHeight);
    const QString manMsg    = QString("   Tab Before Lock x=%1 y=%2 w=%3 h=%4").arg(m_tabViewportManualLayoutBeforeAspectLocking[0]).arg(m_tabViewportManualLayoutBeforeAspectLocking[1]).arg(m_tabViewportManualLayoutBeforeAspectLocking[2]).arg(m_tabViewportManualLayoutBeforeAspectLocking[3]);
    AString msgOut;
    if (m_chartDataViewportValidFlag) {
        const QString chartProjectionMsg = QString("   Chart Projection=" + m_chartDataProjectionMatrix.toFormattedString("      "));
        const QString chartModelMsg = QString("   Chart Model View=" + m_chartDataModelViewMatrix.toFormattedString("      "));
        const QString chartViewportMsg = QString("   Chart x=%1 y=%2 w=%3 h=%4").arg(m_chartDataX).arg(m_chartDataY).arg(m_chartDataWidth).arg(m_chartDataHeight);
        msgOut.appendWithNewLine(chartProjectionMsg);
        msgOut.appendWithNewLine(chartModelMsg);
        msgOut.appendWithNewLine(chartViewportMsg);
    }
    else {
        msgOut.appendWithNewLine("   Chart Invalid.");
    }
    msgOut.appendWithNewLine(windowBeforeLockMsg);
    msgOut.appendWithNewLine(windowMsg);
    msgOut.appendWithNewLine(tabMsg);
    msgOut.appendWithNewLine(modelMsg);
    msgOut.appendWithNewLine(manMsg);
    
    return msgOut;
}

/**
 * Create viewport contents for a single tab using the given window content and window sizes.
 *
 * @param browserTabContent
 *     Tab's content that is being drawn.
 * @param gapsAndMargins
 *     Gaps and margins
 * @param windowIndex
 *    Index of browser window.
 * @param windowBeforeAspectLockingViewport
 *    Viewport of window before aspect locking is applied
 * @param windowViewport
 *    Viewport of WINDOW in which drawing takes place.
 * @return
 *     Viewport content for the single tab.
 */
BrainOpenGLViewportContent*
BrainOpenGLViewportContent::createViewportForSingleTab(std::vector<BrowserTabContent*>& allTabContents,
                                                       BrowserTabContent* selectedTabContent,
                                                       const GapsAndMargins* gapsAndMargins,
                                                       const int32_t windowIndex,
                                                       const int32_t windowBeforeAspectLockingViewport[4],
                                                       const int32_t windowViewport[4])
{
    const int tabViewportBeforeAspectLocking[4] = {
        windowViewport[0],
        windowViewport[1],
        windowViewport[2],
        windowViewport[3]
    };
    
    int tabViewport[4] = {
        windowViewport[0],
        windowViewport[1],
        windowViewport[2],
        windowViewport[3]
    };
    
    int modelViewport[4] = {
        tabViewport[0],
        tabViewport[1],
        tabViewport[2],
        tabViewport[3]
    };
    
    
    std::unique_ptr<EventBrowserWindowContent> eventContent = EventBrowserWindowContent::getWindowContent(windowIndex);
    EventManager::get()->sendEvent(eventContent->getPointer());
    
    const BrowserWindowContent* browserWindowContent = eventContent->getBrowserWindowContent();
    CaretAssert(browserWindowContent);
    if (browserWindowContent->isAllTabsInWindowAspectRatioLocked()) {
        /**
         * Update aspect locking for any tabs that have not bee locked
         * Locking is done here so that it will work for new tabs and
         * for tabs from old scenes before "lock all".
         * Aspect locking is also performed here so that it works 
         * with both GUI and Command Line Show Scene
         */
        for (auto btc : allTabContents) {
            if ( ! btc->isAspectRatioLocked()) {
                if (tabViewport[2] > 0) {
                    const float aspectRatio = (static_cast<float>(tabViewport[3])
                                   / static_cast<float>(tabViewport[2]));
                    
                    btc->setAspectRatio(aspectRatio);
                    btc->setAspectRatioLocked(true);
                }
            }
        }
    }
    
    if (selectedTabContent != NULL) {
        if (selectedTabContent->isAspectRatioLocked()) {
            const float aspectRatio = selectedTabContent->getAspectRatio();
            BrainOpenGLViewportContent::adjustViewportForAspectRatio(tabViewport,
                                                                     aspectRatio);
        }
        
        const int32_t tabIndex = selectedTabContent->getTabNumber();
        createModelViewport(tabViewport,
                            tabIndex,
                            gapsAndMargins,
                            modelViewport);
    }
    

    
    BrainOpenGLViewportContent* vpContent = new BrainOpenGLViewportContent(windowBeforeAspectLockingViewport,
                                                                           windowViewport,
                                                                           tabViewportBeforeAspectLocking,
                                                                           tabViewport,
                                                                           modelViewport,
                                                                           windowIndex,
                                                                           false,
                                                                           selectedTabContent,
                                                                           NULL); /* no spacer tab */
    
    return vpContent;
}

/**
 * Create Viewport Contents for the given tab contents, window sizes, and tile sizes.
 *
 * @param tabContents
 *     Content of each tab.
 * @param browserWindowContent
 *     Content of window.
 * @param gapsAndMargins
 *     Contains margins around edges of tabs
 * @param windowBeforeAspectLockingViewport
 *    Viewport of window before aspect locking is applied
 * @param windowViewport
 *     The window's viewport.
 * @param windowIndex
 *     Index of the window.
 * @param hightlightTabIndex
 *     Index of tab that is highlighted when selected by user.
 * @return
 *     Vector containing data for drawing each model.
 */
std::vector<BrainOpenGLViewportContent*>
BrainOpenGLViewportContent::createViewportContentForTileTabs(std::vector<BrowserTabContent*>& tabContents,
                                                             BrowserWindowContent* browserWindowContent,
                                                             const GapsAndMargins* gapsAndMargins,
                                                             const int32_t windowBeforeAspectLockingViewport[4],
                                                             const int32_t windowViewport[4],
                                                             const int32_t windowIndex,
                                                             const int32_t highlightTabIndex)
{
    CaretAssert(browserWindowContent);
    CaretAssert(gapsAndMargins);
    
    std::vector<BrainOpenGLViewportContent*> viewportContentsOut;
    
    const int32_t numberOfTabs = static_cast<int32_t>(tabContents.size());
    if (numberOfTabs <= 0) {
        return viewportContentsOut;
    }
    
    switch (browserWindowContent->getTileTabsConfigurationMode()) {
        case TileTabsLayoutConfigurationTypeEnum::AUTOMATIC_GRID:
        case TileTabsLayoutConfigurationTypeEnum::CUSTOM_GRID:
        {
            TileTabsLayoutBaseConfiguration* tileTabsConfiguration = browserWindowContent->getSelectedTileTabsGridConfiguration();
            CaretAssert(tileTabsConfiguration);
            viewportContentsOut = createViewportContentForGridTileTabs(tabContents,
                                                                       browserWindowContent,
                                                                       tileTabsConfiguration->castToGridConfiguration(),
                                                                       gapsAndMargins,
                                                                       windowBeforeAspectLockingViewport,
                                                                       windowViewport,
                                                                       windowIndex,
                                                                       highlightTabIndex);
        }
            break;
        case TileTabsLayoutConfigurationTypeEnum::MANUAL:
            viewportContentsOut = createViewportContentForManualTileTabs(tabContents,
                                                                         browserWindowContent,
                                                                         gapsAndMargins,
                                                                         windowBeforeAspectLockingViewport,
                                                                         windowViewport,
                                                                         windowIndex,
                                                                         highlightTabIndex);
            break;
    }

    return viewportContentsOut;
}

/**
 * Create Viewport Contents for the given tab contents, window sizes, and tile sizes
 * for a grid configuration.
 *
 * @param tabContents
 *     Content of each tab.
 * @param browserWindowContent
 *     Content of window.
 * @param gridConfiguration
 *     The grid configuration
 * @param gapsAndMargins
 *     Contains margins around edges of tabs
 * @param windowBeforeAspectLockingViewport
 *    Viewport of window before aspect locking is applied
 * @param windowViewport
 *     The window's viewport.
 * @param windowIndex
 *     Index of the window.
 * @param hightlightTabIndex
 *     Index of tab that is highlighted when selected by user.
 * @return
 *     Vector containing data for drawing each model.
 */
std::vector<BrainOpenGLViewportContent*>
BrainOpenGLViewportContent::createViewportContentForGridTileTabs(std::vector<BrowserTabContent*>& tabContents,
                                                                 BrowserWindowContent* browserWindowContent,
                                                                 TileTabsLayoutGridConfiguration* gridConfiguration,
                                                                 const GapsAndMargins* gapsAndMargins,
                                                                 const int32_t windowBeforeAspectLockingViewport[4],
                                                                 const int32_t windowViewport[4],
                                                                 const int32_t windowIndex,
                                                                 const int32_t highlightTabIndex)
{
    std::vector<BrainOpenGLViewportContent*> viewportContentsOut;
    
    CaretAssert(browserWindowContent);
    CaretAssert(gapsAndMargins);
    CaretAssert(gridConfiguration);
    
    int32_t windowX      = windowViewport[0];
    int32_t windowY      = windowViewport[1];
    const int32_t windowWidth  = windowViewport[2];
    const int32_t windowHeight = windowViewport[3];
    
    /*
     * The tile tabs configuration provides the sizes of the 
     * rows and columns.  The user may "stretch" rows and/or
     * columns.  Thus, the tabs viewports may not be uniformly sized.
     */
    std::vector<int32_t> rowHeights;
    std::vector<int32_t> columnWidths;
    
    const int32_t numberOfTabs = static_cast<int32_t>(tabContents.size());
    gridConfiguration->getRowHeightsAndColumnWidthsForWindowSize(windowWidth,
                                                                     windowHeight,
                                                                     numberOfTabs,
                                                                     browserWindowContent->getTileTabsConfigurationMode(),
                                                                     rowHeights,
                                                                     columnWidths);
    
    const int32_t numRows    = static_cast<int32_t>(rowHeights.size());
    const int32_t numColumns = static_cast<int32_t>(columnWidths.size());
    const int32_t numCells   = numRows * numColumns;
    if (numCells <= 0) {
        return viewportContentsOut;
    }
    CaretAssert(numRows > 0);
    CaretAssert(numColumns > 0);
    
    const bool allTabsAspectLockedFlag = browserWindowContent->isAllTabsInWindowAspectRatioLocked();
    
    /*
     * Create the sizes for the tabs before and after application
     * of aspect ratio.
     */
    std::vector<TileTabsViewportSizingInfo> tabSizeInfoVector;
    int32_t iTab = 0;
    for (int32_t iRowFromTop = 0; iRowFromTop < numRows; iRowFromTop++) {
        CaretAssertVectorIndex(rowHeights, iRowFromTop);
        const int32_t vpHeight = rowHeights[iRowFromTop];

        for (int32_t jCol = 0; jCol < numColumns; jCol++) {
            bool spacerTabFlag = false;
            const TileTabsGridRowColumnContentTypeEnum::Enum rowContentType = gridConfiguration->getRow(iRowFromTop)->getContentType();
            switch (rowContentType) {
                case TileTabsGridRowColumnContentTypeEnum::SPACE:
                    spacerTabFlag = true;
                    break;
                case TileTabsGridRowColumnContentTypeEnum::TAB:
                    break;
            }
            
            const TileTabsGridRowColumnContentTypeEnum::Enum tabContentType = gridConfiguration->getColumn(jCol)->getContentType();
            switch (tabContentType) {
                case TileTabsGridRowColumnContentTypeEnum::SPACE:
                    spacerTabFlag = true;
                    break;
                case TileTabsGridRowColumnContentTypeEnum::TAB:
                    break;
            }

            CaretAssertVectorIndex(columnWidths, jCol);
            const int32_t vpWidth = columnWidths[jCol];
            
            if (spacerTabFlag) {
                EventSpacerTabGet spacerTabEvent(windowIndex, iRowFromTop, jCol);
                EventManager::get()->sendEvent(spacerTabEvent.getPointer());
                SpacerTabContent* spacerTabContent = spacerTabEvent.getSpacerTabContent();
                if (spacerTabContent != NULL) {
                    BrowserTabContent* invalidBrowserTab(NULL);
                    TileTabsViewportSizingInfo tsi(invalidBrowserTab,
                                                   spacerTabContent,
                                                   iRowFromTop,
                                                   jCol,
                                                   vpWidth,
                                                   vpHeight);
                    tabSizeInfoVector.push_back(tsi);
                }
                else {
                    AString msg("Failed to get SpacerTabContent for windowIndex=%1, row=%2, column=%3");
                    msg = msg.arg(windowIndex).arg(iRowFromTop).arg(jCol);
                    CaretLogSevere(msg);
                }
            }
            else if (iTab < numberOfTabs) {
                CaretAssertVectorIndex(tabContents,
                                       iTab);
                
                /*
                 * Is aspect ratio locked for all tabs ?
                 */
                if (allTabsAspectLockedFlag) {
                    /**
                     * Update aspect locking for any tabs that have not been locked
                     * Locking is done here so that it will work for new tabs and
                     * for tabs from old scenes before "lock all"
                     * Aspect locking is also performed here so that it works
                     * with both GUI and Command Line Show Scene
                     */
                    if ( ! tabContents[iTab]->isAspectRatioLocked()) {
                        if (vpWidth > 0) {
                            const float aspectRatio = (static_cast<float>(vpHeight)
                                                       / static_cast<float>(vpWidth));
                            
                            tabContents[iTab]->setAspectRatio(aspectRatio);
                            tabContents[iTab]->setAspectRatioLocked(true);
                        }
                    }
                }
                
                /*
                 * Note: the constructor will adjust the width and
                 * height if lock aspect ratio is enabled.
                 */
                SpacerTabContent* invalidSpacerTab(NULL);
                TileTabsViewportSizingInfo tsi(tabContents[iTab],
                                               invalidSpacerTab,
                                               iRowFromTop,
                                               jCol,
                                               vpWidth,
                                               vpHeight);
                
                tabSizeInfoVector.push_back(tsi);
                iTab++;
            }
        }
        
    }
    
    if (gridConfiguration->isCenteringCorrectionEnabled()) {
        /*
         * Kludge that fixes old scenes
         * NEED TO ADJUST FOR X too ?
         * THIS WILL NEED TO BE AN OPTION SOMEWHERE
         * TRY MORE THAN TWO ROWS
         */
        std::vector<int32_t> tabRowHeights(numRows, 0);
        std::vector<int32_t> tabColumnWidths(numColumns, 0);
        {
            for (auto tsi : tabSizeInfoVector) {
                const int32_t rowIndex = tsi.m_rowIndexFromTop;
                tabRowHeights[rowIndex] = std::max(tabRowHeights[rowIndex],
                                                   tsi.m_height);
                
                const int32_t columnIndex(tsi.m_columnIndex);
                tabColumnWidths[columnIndex] = std::max(tabColumnWidths[columnIndex],
                                                        tsi.m_width);
            }
            
            for (int32_t i = 0; i < numRows; i++) {
                rowHeights[i] = tabRowHeights[i];
            }
            for (int32_t i = 0; i < numColumns; i++) {
                columnWidths[i] = tabColumnWidths[i];
            }
        }
        const int32_t allTabsHeight = std::accumulate(tabRowHeights.begin(), tabRowHeights.end(), 0);
        windowY -= ((windowHeight - allTabsHeight) / 2);
        
        const int32_t allTabsWidth = std::accumulate(tabColumnWidths.begin(), tabColumnWidths.end(), 0);
        const int32_t offset = (windowWidth - allTabsWidth) / 2;
        windowX += offset;
    }
    
    
    
    /*
     * Note: There may be more tabs than there are cells (rows * columns)
     * so some tabs may not be displayed.
     */
    const int32_t numberOfDisplayedTabs = static_cast<int32_t>(tabSizeInfoVector.size());
    
    /*
     * Set the X and Y-coordinates for the tab viewports
     * We start at the bottom row, left corner.
     */
    int32_t vpY = windowHeight + windowY;
    for (int32_t iRow = 0; iRow < numRows; iRow++) {
        CaretAssertVectorIndex(rowHeights, iRow);
        vpY -= rowHeights[iRow];
        int32_t vpX = windowX;
        for (int32_t jCol = 0; jCol < numColumns; jCol++) {
            
            TileTabsViewportSizingInfo* tabSizePtr = NULL;
            for (int32_t iTab = 0; iTab < numberOfDisplayedTabs; iTab++) {
                CaretAssertVectorIndex(tabSizeInfoVector, iTab);
                if ((tabSizeInfoVector[iTab].m_rowIndexFromTop == iRow)
                    && (tabSizeInfoVector[iTab].m_columnIndex == jCol)) {
                    tabSizePtr = &tabSizeInfoVector[iTab];
                    break;
                }
            }
            
            if (tabSizePtr != NULL) {
                CaretAssertVectorIndex(columnWidths, jCol);
                CaretAssertVectorIndex(rowHeights, iRow);
                int32_t tabX = vpX;
                int32_t tabY = vpY;
                
                /*
                 * Adjust Tab's X-coord so centered in tab
                 * when lock aspect causes tab's width to be smaller
                 * than the column's width
                 */
                CaretAssertVectorIndex(columnWidths, jCol);
                const int32_t extraWidth = columnWidths[jCol] - tabSizePtr->m_width;
                if (extraWidth > 1) {
                    const int32_t halfExtraWidth = extraWidth / 2;
                    tabX += halfExtraWidth;
                }
                
                /*
                 * Adjust Tab's Y-coord so centered in tab
                 * when lock aspect causes tab's height to be smaller
                 * than the row's height
                 */
                CaretAssertVectorIndex(rowHeights, iRow);
                const int32_t extraHeight = rowHeights[iRow] - tabSizePtr->m_height;
                if (extraHeight > 1) {
                    const int32_t halfExtraHeight = extraHeight / 2;
                    tabY += halfExtraHeight;
                }
                
                const int tabViewport[4] = {
                    tabX,
                    tabY,
                    tabSizePtr->m_width,
                    tabSizePtr->m_height
                };

                /*
                 * Model is drawn in the model viewport inside any margins.
                 */
                int32_t tabIndex = -1;
                bool highlightTabFlag = false;
                if (tabSizePtr->m_browserTabContent != NULL) {
                    tabIndex = tabSizePtr->m_browserTabContent->getTabNumber();
                    highlightTabFlag = (highlightTabIndex == tabIndex);
                }
                int modelViewport[4] = { 0, 0, 0, 0 };
                createModelViewport(tabViewport,
                                    tabIndex,
                                    gapsAndMargins,
                                    modelViewport);

                BrainOpenGLViewportContent* vpContent = new BrainOpenGLViewportContent(windowBeforeAspectLockingViewport,
                                                                                       windowViewport,
                                                                                       tabViewport,
                                                                                       tabViewport,
                                                                                       modelViewport,
                                                                                       browserWindowContent->getWindowIndex(),
                                                                                       highlightTabFlag,
                                                                                       tabSizePtr->m_browserTabContent,
                                                                                       tabSizePtr->m_spacerTabContent);
                viewportContentsOut.push_back(vpContent);
            }
            else {
                /*
                 * If the number of tabs is less than the number of tiles,
                 * empty viewport content is needed so that user can draw
                 * annotations in these regions.
                 */
                const int tabViewport[4] = {
                    vpX,
                    vpY,
                    columnWidths[jCol],
                    rowHeights[iRow]
                };
                BrainOpenGLViewportContent* vpContent = new BrainOpenGLViewportContent(windowBeforeAspectLockingViewport,
                                                                                       windowViewport,
                                                                                       tabViewport,
                                                                                       tabViewport,
                                                                                       tabViewport,
                                                                                       browserWindowContent->getWindowIndex(),
                                                                                       false,
                                                                                       NULL,
                                                                                       NULL);
                viewportContentsOut.push_back(vpContent);
            }
            
            CaretAssertVectorIndex(columnWidths, jCol);
            vpX += columnWidths[jCol];
        }
    }
 
    const bool printViewportsFlag(false);
    if (printViewportsFlag) {
        for (auto vpc : viewportContentsOut) {
            if (vpc->getTabIndex() == 0) {
                std::cout << "TAB 1: " << vpc->toString() << std::endl << std::endl;
            }
        }
    }
    
    return viewportContentsOut;
}

/**
 * Create Viewport Contents for the given tab contents, window sizes, and tile sizes
 * for a manual configuration.
 *
 * @param tabContents
 *     Content of each tab.
 * @param browserWindowContent
 *     Content of window.
 * @param gapsAndMargins
 *     Contains margins around edges of tabs
 * @param windowBeforeAspectLockingViewport
 *    Viewport of window before aspect locking is applied
 * @param windowViewport
 *     The window's viewport.
 * @param windowIndex
 *     Index of the window.
 * @param hightlightTabIndex
 *     Index of tab that is highlighted when selected by user.
 * @return
 *     Vector containing data for drawing each model.
 */
std::vector<BrainOpenGLViewportContent*>
BrainOpenGLViewportContent::createViewportContentForManualTileTabs(std::vector<BrowserTabContent*>& tabContents,
                                                                   BrowserWindowContent* browserWindowContent,
                                                                   const GapsAndMargins* gapsAndMargins,
                                                                   const int32_t windowBeforeAspectLockingViewport[4],
                                                                   const int32_t windowViewport[4],
                                                                   const int32_t /*windowIndex*/,
                                                                   const int32_t highlightTabIndex)
{
    const bool allTabsAspectLockedFlag = browserWindowContent->isAllTabsInWindowAspectRatioLocked();

    std::vector<BrainOpenGLViewportContent*> viewportContentsOut;
    
    CaretAssert(browserWindowContent);
    CaretAssert(gapsAndMargins);
    
    int32_t windowX      = windowViewport[0];
    int32_t windowY      = windowViewport[1];
    const int32_t windowWidth  = windowViewport[2];
    const int32_t windowHeight = windowViewport[3];

    if (tabContents.size() > 1) {
        /*
         * Sort tabs by stacking order
         */
        std::sort(tabContents.begin(),
                  tabContents.end(),
                  [](BrowserTabContent* a, BrowserTabContent* b) { return (a->getManualLayoutBrowserTabAnnotation()->getStackingOrder()
                                                                           > b->getManualLayoutBrowserTabAnnotation()->getStackingOrder()); } );
    }
    
    for (const auto tab : tabContents) {
        CaretAssert(tab);
        
        const AnnotationBrowserTab* tabAnnotation = tab->getManualLayoutBrowserTabAnnotation();
        if ( ! tabAnnotation->isBrowserTabDisplayed()) {
            continue;
        }
        
        /*
         * Tab Geometry is in percentages ranging [0.0, 100.0]
         */
        float minX(0.0), maxX(0.0), minY(0.0), maxY(0.0);
        tabAnnotation->getBounds2D(minX, maxX, minY, maxY);
        const float tabAnnWidth(maxX - minX);
        const float tabAnnHeight(maxY - minY);
        int32_t tabX(((minX / 100.0) * windowWidth) + windowX);
        int32_t tabY(((minY / 100.0) * windowHeight) + windowY);
        int32_t tabWidth((tabAnnWidth / 100.0) * windowWidth);
        int32_t tabHeight((tabAnnHeight / 100.0) * windowHeight);

        /*
         * Is aspect ratio locked for all tabs ?
         */
        if (allTabsAspectLockedFlag) {
            /**
             * Update aspect locking for any tabs that have not been locked
             * Locking is done here so that it will work for new tabs and
             * for tabs from old scenes before "lock all"
             * Aspect locking is also performed here so that it works
             * with both GUI and Command Line Show Scene
             */
            if ( ! tab->isAspectRatioLocked()) {
                if (tabWidth > 0) {
                    const float aspectRatio = (static_cast<float>(tabHeight)
                                               / static_cast<float>(tabWidth));
                    
                    tab->setAspectRatio(aspectRatio);
                    tab->setAspectRatioLocked(true);
                }
            }
        }
        
        int tabViewportBeforeAspectLocking[4] {
            tabX,
            tabY,
            tabWidth,
            tabHeight
        };
        
        /*
         * Adjust size of tab to match aspect ratio
         */
        if (tab->isAspectRatioLocked()) {
            const int32_t nonLockedWidth(tabWidth);
            const int32_t nonLockedHeight(tabHeight);
            BrainOpenGLViewportContent::adjustWidthHeightForAspectRatio(tab->getAspectRatio(),
                                                                        tabWidth,
                                                                        tabHeight);
            
            /*
             * Center tab within its region
             */
            const int32_t offsetX = nonLockedWidth - tabWidth;
            if (offsetX > 0) {
                tabX += (offsetX / 2);
            }
            
            const int32_t offsetY = nonLockedHeight - tabHeight;
            if (offsetY > 0) {
                tabY += (offsetY / 2);
            }
        }
        
        /*
         * Tab viewport may be invalid (width or height) when the user
         * is editing the layout in tile tags configuration dialog
         */
        if ((tabWidth > 0)
            && (tabHeight > 0)) {
            const int tabViewport[4] = {
                tabX,
                tabY,
                tabWidth,
                tabHeight
            };
            
            /*
             * Model is drawn in the model viewport inside any margins.
             */
            const int32_t tabIndex = tab->getTabNumber();
            const bool highlightTabFlag(highlightTabIndex == tabIndex);
            int modelViewport[4] = { 0, 0, 0, 0 };
            createModelViewport(tabViewport,
                                tabIndex,
                                gapsAndMargins,
                                modelViewport);

            SpacerTabContent* invalidSpaceTabContent(NULL);
            BrainOpenGLViewportContent* vpContent = new BrainOpenGLViewportContent(windowBeforeAspectLockingViewport,
                                                                                   windowViewport,
                                                                                   tabViewportBeforeAspectLocking,
                                                                                   tabViewport,
                                                                                   modelViewport,
                                                                                   browserWindowContent->getWindowIndex(),
                                                                                   highlightTabFlag,
                                                                                   tab,
                                                                                   invalidSpaceTabContent);
            viewportContentsOut.push_back(vpContent);
        }
    }
    
    return viewportContentsOut;
}

/**
 * Create viewport from the model using the tab's viewport and the tabs margins.
 *
 * @paramj tabViewport
 *     Viewport for the tab.
 * @param tabIndex
 *     Index of the tab.
 * @param gapsAndMargins
 *     Gaps and margins.
 * @param modelViewportOut
 *     Output containing viewport for model.
 */
void
BrainOpenGLViewportContent::createModelViewport(const int tabViewport[4],
                                                const int32_t tabIndex,
                                                const GapsAndMargins* gapsAndMargins,
                                                int modelViewportOut[4])
{
    int32_t leftMargin   = 0;
    int32_t rightMargin  = 0;
    int32_t bottomMargin = 0;
    int32_t topMargin    = 0;
    
    modelViewportOut[0] = tabViewport[0];
    modelViewportOut[1] = tabViewport[1];
    modelViewportOut[2] = tabViewport[2];
    modelViewportOut[3] = tabViewport[3];
    
    if (tabIndex >= 0) {
        if (gapsAndMargins != NULL) {
            gapsAndMargins->getMarginsInPixelsForDrawing(tabIndex,
                                                         tabViewport[2],
                                                         tabViewport[3],
                                                         leftMargin,
                                                         rightMargin,
                                                         bottomMargin,
                                                         topMargin);
            const int32_t marginHorizSize = (leftMargin   + rightMargin);
            const int32_t marginVertSize  = (bottomMargin + topMargin);
            if ((marginHorizSize < modelViewportOut[2])
                && (marginVertSize < modelViewportOut[3])) {
                modelViewportOut[0] += leftMargin;
                modelViewportOut[1] += bottomMargin;
                modelViewportOut[2] -= marginHorizSize;
                modelViewportOut[3] -= marginVertSize;
            }
            else {
                CaretLogSevere("Margins are too big for tab "
                               + AString::number(tabIndex + 1)
                               + " viewport.  Viewport (x,y,w,h)="
                               + AString::fromNumbers(modelViewportOut, 4, ",")
                               + " margin (l,r,b,t)="
                               + AString::number(leftMargin) + ","
                               + AString::number(rightMargin) + ","
                               + AString::number(bottomMargin) + ","
                               + AString::number(topMargin));
            }
        }
    }
}

/**
 * When viewing a surface montage, "getModelViewport()" returns the viewport
 * for the entire surface montage (viewport containing all of the surfaces).
 * However, the are instance where we need the viewport for one of surfaces
 * in a surface montage.  This method will find the viewport for one of the 
 * surfaces within the surface montage using the given X and Y coordinates.
 *
 * @param montageX
 *     X-coordinate within the surface montage.
 * @param montageX
 *     X-coordinate within the surface montage.
 * @param subViewportOut
 *     Viewport for individual surface in the montage at the given X, Y 
 *     coordinates.   Note: if a surface montage is not displayed, 
 *     this method returns the same viewport as getModelViewport().
 */
void
BrainOpenGLViewportContent::getSurfaceMontageModelViewport(const int32_t montageX,
                                                           const int32_t montageY,
                                                           int subViewportOut[4]) const
{
    getModelViewport(subViewportOut);
    
    if (this->m_browserTabContent == NULL) {
        return;
    }
    
    ModelSurfaceMontage* msm = m_browserTabContent->getDisplayedSurfaceMontageModel();
    if (msm != NULL) {
        std::vector<const SurfaceMontageViewport*> montageViewports;
        msm->getSurfaceMontageViewportsForTransformation(m_browserTabContent->getTabNumber(),
                                                         montageViewports);
        
        const int x = montageX;
        const int y = montageY;
        
        for (std::vector<const SurfaceMontageViewport*>::const_iterator iter = montageViewports.begin();
             iter != montageViewports.end();
             iter++) {
            const SurfaceMontageViewport* smv = *iter;
            
            int vp[4];
            smv->getViewport(vp);
            
            const int offsetX = x - vp[0];
            if ((offsetX >= 0)
                && (offsetX < vp[2])) {
                const int offsetY = y - vp[1];
                
                if ((offsetY >= 0)
                    && (offsetY < vp[3])) {
                    smv->getViewport(subViewportOut);
                    return;
                }
            }
        }
    }
}

/**
 * Set the object to window transformation for histology.  This instance will take ownership of the given transform and delete when appropriate.
 */
void
BrainOpenGLViewportContent::setHistologyGraphicsObjectToWindowTransform(GraphicsObjectToWindowTransform* transform) const
{
    m_histologyGraphicsObjectToWindowTransform.reset(transform);
}

/**
 * @return Object to window transformation (MAY BE NULL !!!) for histology
 */
const GraphicsObjectToWindowTransform*
BrainOpenGLViewportContent::getHistologyGraphicsObjectToWindowTransform() const
{
    return m_histologyGraphicsObjectToWindowTransform.get();
}

/**
 * Set the object to window transformation for media.  This instance will take ownership of the given transform and delete when appropriate.
 */
void
BrainOpenGLViewportContent::setMediaGraphicsObjectToWindowTransform(GraphicsObjectToWindowTransform* transform) const
{
    m_mediaGraphicsObjectToWindowTransform.reset(transform);
}

/**
 * @return Object to window transformation (MAY BE NULL !!!) for media
 */
const GraphicsObjectToWindowTransform*
BrainOpenGLViewportContent::getMediaGraphicsObjectToWindowTransform() const
{
    return m_mediaGraphicsObjectToWindowTransform.get();
}

/**
 * Set the object to window transformation for volume.  This instance will take ownership of the given transform and delete when appropriate.
 */
void
BrainOpenGLViewportContent::setVolumeMprGraphicsObjectToWindowTransform(const VolumeSliceViewPlaneEnum::Enum sliceViewPlane,
                                                                        GraphicsObjectToWindowTransform* transform) const
{
    switch (sliceViewPlane) {
        case VolumeSliceViewPlaneEnum::ALL:
            CaretAssert(0);
            break;
        case VolumeSliceViewPlaneEnum::AXIAL:
            m_volumeAxialGraphicsObjectToWindowTransform.reset(transform);
            break;
        case VolumeSliceViewPlaneEnum::CORONAL:
            m_volumeCoronalGraphicsObjectToWindowTransform.reset(transform);
            break;
        case VolumeSliceViewPlaneEnum::PARASAGITTAL:
            m_volumeParasagittalGraphicsObjectToWindowTransform.reset(transform);
            break;
    }
}

/**
 * @return Object to window transformation (MAY BE NULL !!!) for volume
 */
const GraphicsObjectToWindowTransform*
BrainOpenGLViewportContent::getVolumeGraphicsObjectToWindowTransform(const VolumeSliceViewPlaneEnum::Enum sliceViewPlane) const
{
    switch (sliceViewPlane) {
        case VolumeSliceViewPlaneEnum::ALL:
            CaretAssert(0);
            break;
        case VolumeSliceViewPlaneEnum::AXIAL:
            return m_volumeAxialGraphicsObjectToWindowTransform.get();
            break;
        case VolumeSliceViewPlaneEnum::CORONAL:
            return m_volumeCoronalGraphicsObjectToWindowTransform.get();
            break;
        case VolumeSliceViewPlaneEnum::PARASAGITTAL:
            return m_volumeParasagittalGraphicsObjectToWindowTransform.get();
            break;
    }
    return NULL;
}

/**
 * @return Step value for custom view dialog that depends upon model in selected tab
 */
float
BrainOpenGLViewportContent::getTranslationStepValueForCustomViewDialog() const
{
    float customViewStepValue(1.0);
    float mousePanningFactor(1.0);
    getTranslationFactors(customViewStepValue,
                          mousePanningFactor);
    return customViewStepValue;
}

/**
 * @return Factor for panning with mouse depends upon model in selected tab
 */
float
BrainOpenGLViewportContent::getTranslationFactorForMousePanning() const
{
    float customViewStepValue(1.0);
    float mousePanningFactor(1.0);
    getTranslationFactors(customViewStepValue,
                          mousePanningFactor);
    return mousePanningFactor;
}

/**
 * Get factors for translation
 * @param customViewStepValueOut
 *    Step value for translation on custom view dialog
 * @param mousePanningFactor
 *    Factor for panning with mouse
 */
void
BrainOpenGLViewportContent::getTranslationFactors(float& customViewStepValueOut,
                                                  float& mousePanningFactorOut) const
{
    customViewStepValueOut = 1.0;
    mousePanningFactorOut  = 1.0;
    
    BrowserTabContent* tabContent = getBrowserTabContent();
    if (tabContent->isMediaDisplayed()) {
        switch (tabContent->getMediaDisplayCoordinateMode()) {
            case MediaDisplayCoordinateModeEnum::PIXEL:
            {
                MediaOverlaySet* overlaySet = tabContent->getMediaOverlaySet();
                CaretAssert(overlaySet);
                MediaOverlay* underlay = overlaySet->getBottomMostEnabledOverlay();
                if (underlay != NULL) {
                    const MediaOverlay::SelectionData selectionData(underlay->getSelectionData());
                    if (selectionData.m_selectedMediaFile != NULL) {
                        const float height(selectionData.m_selectedMediaFile->getHeight());
                        if (height > 0.0) {
                            int32_t viewport[4];
                            getModelViewport(viewport);
                            const float viewportHeight(viewport[3]);
                            if (viewportHeight > 0.0) {
                                const float factor(height / viewportHeight);
                                mousePanningFactorOut = factor * 0.5;
                                customViewStepValueOut = factor;
                            }
                        }
                    }
                }
            }
                break;
            case MediaDisplayCoordinateModeEnum::PLANE:
                break;
        }
    }
}





/* =================================================================================================== */


/**
 * Contructor for tile tabs viewport sizing.
 *
 * @param browserTabContent
 *    Content of the browser tab.
 * @param rowIndexFromTop
 *    Row index starting with top row
 * @param columnIndex
 *    Column index starting on left.
 * @param initialWidth
 *    Initial width of the tab prior to application of aspect ratio.
 * @param initialHeight
 *    Initial height of the tab prior to application of aspect ratio.
 */
BrainOpenGLViewportContent::TileTabsViewportSizingInfo::TileTabsViewportSizingInfo(BrowserTabContent* browserTabContent,
                                                                                   SpacerTabContent* spacerTabContent,
                                                                                   const int32_t rowIndexFromTop,
                                                                                   const int32_t columnIndex,
                                                                                   const float initialWidth,
                                                                                   const float initialHeight)
: m_browserTabContent(browserTabContent),
m_spacerTabContent(spacerTabContent),
m_rowIndexFromTop(rowIndexFromTop),
m_columnIndex(columnIndex),
m_initialWidth(initialWidth),
m_initialHeight(initialHeight),
m_width(initialWidth),
m_height(initialHeight)
{
    if (m_browserTabContent != NULL) {
        if (m_browserTabContent->isAspectRatioLocked()) {
            const float aspectRatio = m_browserTabContent->getAspectRatio();
            if (aspectRatio > 0.0) {
                BrainOpenGLViewportContent::adjustWidthHeightForAspectRatio(aspectRatio,
                                                                            m_width,
                                                                            m_height);
            }
        }
    }
}

/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
BrainOpenGLViewportContent::TileTabsViewportSizingInfo&
BrainOpenGLViewportContent::TileTabsViewportSizingInfo::operator=(const TileTabsViewportSizingInfo& obj)
{
    if (this != &obj) {
        m_browserTabContent = obj.m_browserTabContent;
        m_spacerTabContent  = obj.m_spacerTabContent;
        m_rowIndexFromTop   = obj.m_rowIndexFromTop;
        m_columnIndex       = obj.m_columnIndex;
        m_initialWidth      = obj.m_initialWidth;
        m_initialHeight     = obj.m_initialHeight;
        m_width             = obj.m_initialWidth;
        m_height            = obj.m_initialHeight;
    }
    
    return *this;
}

/**
 * Print for debugging.
 *
 * @param x
 *     X-coordinate of tab viewport.
 * @param y
 *     Y-coordinate of tab viewport.
 */
void
BrainOpenGLViewportContent::TileTabsViewportSizingInfo::print(const int32_t x,
                                                              const int32_t y)
{
    AString name;
    if (m_browserTabContent != NULL) {
        name = m_browserTabContent->getTabName();
    }
    else if (m_spacerTabContent != NULL) {
        name = m_spacerTabContent->getTabName();
    }
    const QString msg("Model: " + name
                      + "\n   row/col: " + QString::number(m_rowIndexFromTop) + ", " + QString::number(m_columnIndex)
                      + "\n   x/y: " + QString::number(x) + ", " + QString::number(y)
                      + "\n   width/height: " + QString::number(m_width) + ", " + QString::number(m_height));
    std::cout << qPrintable(msg) << std::endl;
}


/**
 * Get the viewport for a slice in all slices view.
 *
 * @param tabViewport
 *    The viewport for the tab containing all slices.
 * @param sliceViewPlane
 *    The plane for slice drawing.  Note: "ALL" is used for orientation axes in oblique view.
 * @param allPlanesLayout
 *    The layout in ALL slices view.
 * @param viewportOut
 *    Output viewport (region of graphics area) for drawing slices.
 */
void
BrainOpenGLViewportContent::getSliceAllViewViewport(const int32_t tabViewport[4],
                                                       const VolumeSliceViewPlaneEnum::Enum sliceViewPlane,
                                                       const VolumeSliceViewAllPlanesLayoutEnum::Enum allPlanesLayout,
                                                       int32_t viewportOut[4])
{
    const int32_t gap = 2;
    const int32_t tabViewportX      = tabViewport[0];
    const int32_t tabViewportY      = tabViewport[1];
    const int32_t tabViewportWidth  = tabViewport[2];
    const int32_t tabViewportHeight = tabViewport[3];
    
    switch (allPlanesLayout) {
        case VolumeSliceViewAllPlanesLayoutEnum::COLUMN_LAYOUT:
        {
            const int32_t vpHeight = (tabViewportHeight - (gap * 2)) / 3;
            const int32_t vpOffsetY = vpHeight + gap;
            const int32_t vpWidth  = tabViewportWidth;
            
            switch (sliceViewPlane) {
                case VolumeSliceViewPlaneEnum::ALL:
                    viewportOut[0] = 0;
                    viewportOut[1] = 0;
                    viewportOut[2] = 0;
                    viewportOut[3] = 0;
                    break;
                case VolumeSliceViewPlaneEnum::AXIAL:
                    viewportOut[0] = tabViewportX;
                    viewportOut[1] = tabViewportY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::CORONAL:
                    viewportOut[0] = tabViewportX;
                    viewportOut[1] = tabViewportY + vpOffsetY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                    viewportOut[0] = tabViewportX;
                    viewportOut[1] = tabViewportY + (vpOffsetY * 2);
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
            }
        }
            break;
        case VolumeSliceViewAllPlanesLayoutEnum::GRID_LAYOUT:
        {
            const int32_t vpWidth   = (tabViewportWidth  - gap) / 2;
            const int32_t vpHeight  = (tabViewportHeight - gap) / 2;
            const int32_t vpOffsetX = vpWidth  + gap;
            const int32_t vpOffsetY = vpHeight + gap;
            switch (sliceViewPlane) {
                case VolumeSliceViewPlaneEnum::ALL:
                    viewportOut[0] = tabViewportX;
                    viewportOut[1] = tabViewportY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::AXIAL:
                    viewportOut[0] = tabViewportX + vpOffsetX;
                    viewportOut[1] = tabViewportY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::CORONAL:
                    viewportOut[0] = tabViewportX + vpOffsetX;
                    viewportOut[1] = tabViewportY + vpOffsetY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                    viewportOut[0] = tabViewportX;
                    viewportOut[1] = tabViewportY + vpOffsetY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
            }
        }
            break;
        case VolumeSliceViewAllPlanesLayoutEnum::ROW_LAYOUT:
        {
            const int32_t vpWidth   = (tabViewportWidth - (gap * 2)) / 3;
            const int32_t vpOffsetX = vpWidth + gap;
            const int32_t vpHeight  = tabViewportHeight;
            
            switch (sliceViewPlane) {
                case VolumeSliceViewPlaneEnum::ALL:
                    viewportOut[0] = 0;
                    viewportOut[1] = 0;
                    viewportOut[2] = 0;
                    viewportOut[3] = 0;
                    break;
                case VolumeSliceViewPlaneEnum::AXIAL:
                    viewportOut[0] = tabViewportX + (vpOffsetX * 2);
                    viewportOut[1] = tabViewportY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::CORONAL:
                    viewportOut[0] = tabViewportX + vpOffsetX;
                    viewportOut[1] = tabViewportY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
                case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                    viewportOut[0] = tabViewportX;
                    viewportOut[1] = tabViewportY;
                    viewportOut[2] = vpWidth;
                    viewportOut[3] = vpHeight;
                    break;
            }
        }
            break;
    }
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
 * @param allPlanesLayout
 *    The layout in ALL slices view.
 * @param sliceViewportOut
 *    Output viewport (region of graphics area) for drawing slices.
 */
VolumeSliceViewPlaneEnum::Enum
BrainOpenGLViewportContent::getSliceViewPlaneForVolumeAllSliceView(const GraphicsViewport& viewport,
                                                                   const VolumeSliceViewAllPlanesLayoutEnum::Enum allPlanesLayout,
                                                                   const Vector3D& mousePressXY,
                                                                   GraphicsViewport& sliceViewportOut)
{
    int32_t sliceViewport[4];
    const VolumeSliceViewPlaneEnum::Enum viewPlane(getSliceViewPlaneForVolumeAllSliceView(viewport.getViewport().data(),
                                                                                          allPlanesLayout,
                                                                                          mousePressXY[0],
                                                                                          mousePressXY[1],
                                                                                          sliceViewport));
    sliceViewportOut = GraphicsViewport(sliceViewport);
    return viewPlane;
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
 * @param allPlanesLayout
 *    The layout in ALL slices view.
 * @param sliceViewportOut
 *    Output viewport (region of graphics area) for drawing slices.
 */
VolumeSliceViewPlaneEnum::Enum
BrainOpenGLViewportContent::getSliceViewPlaneForVolumeAllSliceView(const int32_t viewport[4],
                                                                   const VolumeSliceViewAllPlanesLayoutEnum::Enum allPlanesLayout,
                                                                   const int32_t mousePressX,
                                                                   const int32_t mousePressY,
                                                                   int32_t sliceViewportOut[4])
{
    VolumeSliceViewPlaneEnum::Enum view = VolumeSliceViewPlaneEnum::ALL;
    
    std::vector<VolumeSliceViewPlaneEnum::Enum> allSlicePlanes;
    VolumeSliceViewPlaneEnum::getAllEnums(allSlicePlanes);
    for (auto slicePlane : allSlicePlanes) {
        BrainOpenGLViewportContent::getSliceAllViewViewport(viewport,
                                                           slicePlane,
                                                           allPlanesLayout,
                                                           sliceViewportOut);
        const int32_t vpX = mousePressX - sliceViewportOut[0];
        const int32_t vpY = mousePressY - sliceViewportOut[1];
        if ((vpX >= 0)
            && (vpY >= 0)
            && (vpX < sliceViewportOut[2])
            && (vpY < sliceViewportOut[3])) {
            return slicePlane;
        }
    }
    
    CaretLogSevere("Failed to find slice plane in all sliced view");
    return view;
}

/**
 * @return Pair with GraphicsViewport and slice plane.
 *
 * The viewport at the mouse XY.  If this is a volume slice with ALL view, the
 * viewport will be for the individual slice (axial, coronal, parasagittal) NOT all three slices.
 * Otherwise, the viewport is the same as getModelViewport().
 *
 * The graphics viewport will be invalid if there is a failure to find the viewport.
 */
std::pair<GraphicsViewport,
VolumeSliceViewPlaneEnum::Enum>
BrainOpenGLViewportContent::getVolumeSliceViewportAtMouseXY(const int32_t mouseX,
                                                            const int32_t mouseY) const
{
    std::array<int32_t, 4> viewport;
    getModelViewport(viewport.data());
    
    VolumeSliceViewPlaneEnum::Enum slicePlaneOut = VolumeSliceViewPlaneEnum::ALL;
    const BrowserTabContent* btc(getBrowserTabContent());
    if (btc != NULL) {
        if (btc->isVolumeSlicesDisplayed()) {
            switch (btc->getVolumeSliceViewPlane()) {
                case VolumeSliceViewPlaneEnum::ALL:
                {
                    std::array<int32_t, 4>  sliceViewport;
                    slicePlaneOut = getSliceViewPlaneForVolumeAllSliceView(viewport.data(),
                                                                           btc->getVolumeSlicePlanesAllViewLayout(),
                                                                           mouseX,
                                                                           mouseY,
                                                                           sliceViewport.data());
                    if (slicePlaneOut == VolumeSliceViewPlaneEnum::ALL) {
                        return std::make_pair(GraphicsViewport(),
                                              VolumeSliceViewPlaneEnum::ALL);
                    }
                    viewport = sliceViewport;
                }
                    break;
                case VolumeSliceViewPlaneEnum::AXIAL:
                    break;
                case VolumeSliceViewPlaneEnum::CORONAL:
                    break;
                case VolumeSliceViewPlaneEnum::PARASAGITTAL:
                    break;
            }
        }
    }
    
    return std::make_pair(GraphicsViewport(viewport),
                          slicePlaneOut);
    
}
