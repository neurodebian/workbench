#ifndef __USER_INPUT_MODE_VIEW__H_
#define __USER_INPUT_MODE_VIEW__H_

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
#include "UserInputModeAbstract.h"

namespace caret {

    class BrainOpenGLViewportContent;
    class ChartTwoOverlay;

    class UserInputModeView : public UserInputModeAbstract {
        
    public:
        UserInputModeView(const int32_t windowIndex);
        
        virtual ~UserInputModeView();
        
        virtual void initialize();
        
        virtual void finish();

        virtual void update();
        
        virtual CursorEnum::Enum getCursor() const;
        
        virtual bool keyPressEvent(const KeyEvent& /*keyEvent*/) override;
        
        virtual void mouseLeftDoubleClick(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftClick(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftClickWithShift(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftClickWithCtrlShift(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftDrag(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftDragWithAlt(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftDragWithCtrl(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftDragWithShift(const MouseEvent& mouseEvent);

        virtual void mouseLeftPress(const MouseEvent& mouseEvent);
        
        virtual void mouseLeftRelease(const MouseEvent& mouseEvent);
        
        virtual void gestureEvent(const GestureEvent& gestureEvent);
        
        virtual void showContextMenu(const MouseEvent& mouseEvent,
                                     const QPoint& menuPosition,
                                     BrainOpenGLWidget* openGLWidget);
        
    protected:
        UserInputModeView(const int32_t windowIndex,
                          const UserInputModeEnum::Enum inputMode);
        
    private:
        enum class ChartActiveLayerMode {
            DECREMENT,
            DESELECT_ALL,
            INCREMENT,
            SELECT
        };
        
        UserInputModeView(const UserInputModeView&);

        UserInputModeView& operator=(const UserInputModeView&);

        void updateGraphics(const MouseEvent& mouseEvent);
        
        void updateGraphics(const BrainOpenGLViewportContent* viewportContent);
        
        void processModelViewIdentification(BrainOpenGLViewportContent* viewportContent,
                                            BrainOpenGLWidget* openGLWidget,
                                            const int32_t mouseClickX,
                                            const int32_t mouseClickY);
        
        void processChartActiveLayerAction(const ChartActiveLayerMode chartActiveMode,
                                           ChartTwoOverlay* chartOverlay,
                                           const int32_t pointIndex);
        
        const int32_t m_browserWindowIndex;
        
    public:
        virtual AString toString() const;
        
    };
    
#ifdef __USER_INPUT_MODE_VIEW_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __USER_INPUT_MODE_VIEW_DECLARE__

} // namespace
#endif  //__USER_INPUT_MODE_VIEW__H_
