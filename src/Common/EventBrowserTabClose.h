#ifndef __EVENT_BROWSER_TAB_CLOSE_H__
#define __EVENT_BROWSER_TAB_CLOSE_H__

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

#include <vector>

#include "Event.h"

namespace caret {

    class BrowserTabContent;
    
    /// Event for deleting a browser tab
    class EventBrowserTabClose : public Event {
        
    public:
        EventBrowserTabClose(BrowserTabContent* browserTab,
                             const int32_t browserTabIndex,
                             const int32_t windowIndex,
                             const int32_t windowTabBarPositionIndex);

        virtual ~EventBrowserTabClose();
        
        BrowserTabContent* getBrowserTab();
        
        int32_t getBrowserTabIndex() const;
        
        int32_t getWindowIndex() const;
        
        int32_t getWindowTabBarPositionIndex() const;
        
    private:
        EventBrowserTabClose(const EventBrowserTabClose&);
        
        EventBrowserTabClose& operator=(const EventBrowserTabClose&);
        
        BrowserTabContent* m_browserTab;
        
        const int32_t m_browserTabIndex;
        
        const int32_t m_windowIndex;
        
        const int32_t m_windowTabBarPositionIndex;
    };

} // namespace

#endif // __EVENT_BROWSER_TAB_CLOSE_H__
