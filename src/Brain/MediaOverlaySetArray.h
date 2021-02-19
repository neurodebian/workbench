#ifndef __MEDIA_OVERLAY_SET_ARRAY_H__
#define __MEDIA_OVERLAY_SET_ARRAY_H__

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

#include "EventListenerInterface.h"
#include "Overlay.h"

namespace caret {
    class MediaOverlaySet;
    
    class MediaOverlaySetArray : public CaretObject, public EventListenerInterface {
        
    public:
        MediaOverlaySetArray(const AString& name);
        
        virtual ~MediaOverlaySetArray();
        
        int32_t getNumberOfOverlaySets();
        
        MediaOverlaySet* getOverlaySet(const int32_t indx);
        
        void initializeOverlaySelections();
        
        void copyOverlaySet(const int32_t sourceTabIndex,
                            const int32_t destinationTabIndex);
        
        virtual SceneClass* saveTabIndicesToScene(const std::vector<int32_t>& tabIndices,
                                                  const SceneAttributes* sceneAttributes,
                                                  const AString& instanceName);
        
        virtual void restoreFromScene(const SceneAttributes* sceneAttributes,
                                      const SceneClass* sceneClass);
    private:
        MediaOverlaySetArray(const MediaOverlaySetArray&);

        MediaOverlaySetArray& operator=(const MediaOverlaySetArray&);
        
    public:

        // ADD_NEW_METHODS_HERE

        virtual AString toString() const;
        
        virtual void receiveEvent(Event* event);

    private:
        // ADD_NEW_MEMBERS_HERE
        
        void initialize();

        /** Name for this overlay set array */
        AString m_name;
        
        /** The overlay sets */
        std::vector<MediaOverlaySet*> m_overlaySets;
    };
    
#ifdef __MEDIA_OVERLAY_SET_ARRAY_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __MEDIA_OVERLAY_SET_ARRAY_DECLARE__

} // namespace
#endif  //__MEDIA_OVERLAY_SET_ARRAY_H__
