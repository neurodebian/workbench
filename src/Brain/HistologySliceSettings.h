#ifndef __HISTOLOGY_SLICE_SETTINGS_H__
#define __HISTOLOGY_SLICE_SETTINGS_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2022 Washington University School of Medicine
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



#include <memory>

#include "CaretObject.h"
#include "HistologyCoordinate.h"
#include "SceneableInterface.h"
#include "Vector3D.h"

namespace caret {
    class HistologySlicesFile;
    class SceneClassAssistant;

    class HistologySliceSettings : public CaretObject, public SceneableInterface {
        
    public:
        HistologySliceSettings();
        
        virtual ~HistologySliceSettings();
        
        HistologySliceSettings(const HistologySliceSettings& obj);

        HistologySliceSettings& operator=(const HistologySliceSettings& obj);
        
        void copyYokedSettings(const HistologySlicesFile* histologySlicesFile,
                               const HistologySliceSettings& settings);
        
        void updateForHistologySlicesFile(const HistologySlicesFile* histologySlicesFile);
        
        void reset();
        
        HistologyCoordinate getHistologyCoordinate(const HistologySlicesFile* histologySlicesFile) const;
        
        void setHistologyCoordinate(const HistologyCoordinate& histologyCoordinate);
        
        void selectSlicesAtCenter(const HistologySlicesFile* histologySlicesFile);
        
        // ADD_NEW_METHODS_HERE

        virtual AString toString() const;
        
        virtual SceneClass* saveToScene(const SceneAttributes* sceneAttributes,
                                        const AString& instanceName);

        virtual void restoreFromScene(const SceneAttributes* sceneAttributes,
                                      const SceneClass* sceneClass);

          
          
          
          
          
// If there will be sub-classes of this class that need to save
// and restore data from scenes, these pure virtual methods can
// be uncommented to force their implementation by sub-classes.
//    protected: 
//        virtual void saveSubClassDataToScene(const SceneAttributes* sceneAttributes,
//                                             SceneClass* sceneClass) = 0;
//
//        virtual void restoreSubClassDataFromScene(const SceneAttributes* sceneAttributes,
//                                                  const SceneClass* sceneClass) = 0;

    private:
        void copyHelperHistologySliceSettings(const HistologySliceSettings& obj);

        std::unique_ptr<SceneClassAssistant> m_sceneAssistant;

        mutable HistologyCoordinate m_histologyCoordinate = HistologyCoordinate();
        
        mutable Vector3D m_sliceCoordinateXYZ;
        
        bool m_initializedFlag = false;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __HISTOLOGY_SLICE_SETTINGS_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __HISTOLOGY_SLICE_SETTINGS_DECLARE__

} // namespace
#endif  //__HISTOLOGY_SLICE_SETTINGS_H__
