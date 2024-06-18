#ifndef __MouseEvent_H__
#define __MouseEvent_H__

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

#include "CaretObject.h"
#include "Vector3D.h"

#include <stdint.h>

class QScreen;

namespace caret {
    
    class BrainOpenGLViewportContent;
    class BrainOpenGLWidget;
    class BrainOpenGLWindowContent;
    class BrowserTabContent;
    
    /**
     * Contains information about a mouse event in the OpenGL region.
     */
    class MouseEvent : public CaretObject {
        
    public:
        /**
         * Contains mouse X/Y coordinates
         */
        class XY {
        public:
            XY(const int32_t x,
               const int32_t y)
            : m_x(x),
              m_y(y) { }
            
            int32_t m_x;
            int32_t m_y;
        };
        
        MouseEvent(const BrainOpenGLWindowContent* windowContent,
                   const BrainOpenGLViewportContent* viewportContent,
                   BrainOpenGLWidget* openGLWidget,
                   const int32_t browserWindowIndex,
                   const int32_t x,
                   const int32_t y,
                   const int32_t dx,
                   const int32_t dy,
                   const int32_t mousePressX,
                   const int32_t mousePressY,
                   const std::vector<XY>& mouseHistoryXY,
                   const int32_t globalX,
                   const int32_t globalY,
                   const bool firstDraggingFlag);
        
        virtual ~MouseEvent();
        
        MouseEvent(const MouseEvent& o);
        
        MouseEvent& operator=(const MouseEvent& o);
        
    private:
        void initializeMembersMouseEvent();
        
        void copyHelperMouseEvent(const MouseEvent& me);
        
    public:
        AString toString() const;
        
        BrainOpenGLWindowContent* getWindowContent() const;
        
        BrainOpenGLViewportContent* getViewportContent() const;
        
        BrainOpenGLWidget* getOpenGLWidget() const;

        BrowserTabContent* getBrowserTabContentAsMouseXY() const;
        
        BrowserTabContent* getBrowserTabContentAsMousePressXY() const;
        
        int32_t getBrowserWindowIndex() const;
        
        int32_t getDx() const;
        
        int32_t getDy() const;
        
        int32_t getX() const;
        
        int32_t getY() const;
        
        int32_t getPressedX() const;
        
        int32_t getPressedY() const;
        
        Vector3D getXY() const;
        
        Vector3D getPressedXY() const;
        
        int32_t getXyHistoryCount() const;
        
        XY getHistoryAtIndex(const int32_t index) const;
        
        void offsetXY(const float offsetX,
                      const float offsetY);
        
        int32_t getWheelRotation() const;
        
        QScreen* getScreen() const;
        
        int32_t getGlobalX() const;
        
        int32_t getGlobalY() const;
        
        bool isFirstDragging() const;
        
    private:
        BrainOpenGLWindowContent* m_windowContent;
        
        BrainOpenGLViewportContent* m_viewportContent;
        
        BrainOpenGLWidget* m_openGLWidget;
        
        int32_t m_browserWindowIndex;
        
        int32_t m_x;
        
        int32_t m_y;
        
        int32_t m_dx;
        
        int32_t m_dy;

        int32_t m_pressX;
        
        int32_t m_pressY;
        
        std::vector<XY> m_xyHistory;
        
        int32_t m_wheelRotation;
        
        int32_t m_globalX = -1;
        
        int32_t m_globalY = -1;
        
        bool m_firstDraggingFlag;
    };
    
} // namespace

#endif // __MouseEvent_H__
