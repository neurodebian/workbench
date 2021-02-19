#ifndef __EVENT_GRAPHICS_UPDATE_ONE_WINDOW_H__
#define __EVENT_GRAPHICS_UPDATE_ONE_WINDOW_H__

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


#include "Event.h"

namespace caret {

    /// Event for updating graphics in one window.
    class EventGraphicsUpdateOneWindow : public Event {
        
    public:
        EventGraphicsUpdateOneWindow(const int32_t windowIndex,
                                     const bool doRepaint = false);
        
        virtual ~EventGraphicsUpdateOneWindow();
        
        /// get the index of the window that is to be updated.
        int32_t getWindowIndex() const { return this->windowIndex; }
        
        bool isRepaint() const;
        
    private:
        EventGraphicsUpdateOneWindow(const EventGraphicsUpdateOneWindow&);
        
        EventGraphicsUpdateOneWindow& operator=(const EventGraphicsUpdateOneWindow&);
        
        /** index of window for update */
        int32_t windowIndex;
        
        bool m_doRepaintFlag = false;
    };

} // namespace

#endif // __EVENT_GRAPHICS_UPDATE_ONE_WINDOW_H__
