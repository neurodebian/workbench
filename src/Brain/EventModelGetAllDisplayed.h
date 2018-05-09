#ifndef __EVENT_MODEL_GET_ALL_DISPLAYED_H__
#define __EVENT_MODEL_GET_ALL_DISPLAYED_H__

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

#include <set>

#include "Event.h"


namespace caret {

    class Model;
    
    class EventModelGetAllDisplayed : public Event {
        
    public:
        EventModelGetAllDisplayed();
        
        virtual ~EventModelGetAllDisplayed();
        
        std::set<Model*> getModels() const;
        
        void addModel(Model* model);
        
        // ADD_NEW_METHODS_HERE

    private:
        EventModelGetAllDisplayed(const EventModelGetAllDisplayed&);

        EventModelGetAllDisplayed& operator=(const EventModelGetAllDisplayed&);
        
        std::set<Model*> m_models;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __EVENT_MODEL_GET_ALL_DISPLAYED_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __EVENT_MODEL_GET_ALL_DISPLAYED_DECLARE__

} // namespace
#endif  //__EVENT_MODEL_GET_ALL_DISPLAYED_H__
