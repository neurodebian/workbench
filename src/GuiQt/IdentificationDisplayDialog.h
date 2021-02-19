#ifndef __IDENTIFICATION_DISPLAY_DIALOG_H__
#define __IDENTIFICATION_DISPLAY_DIALOG_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2019 Washington University School of Medicine
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

#include "WuQDialogNonModal.h"

#include "SceneableInterface.h"


namespace caret {
    class IdentificationDisplayWidget;
    class SceneClassAssistant;

    class IdentificationDisplayDialog : public WuQDialogNonModal, public SceneableInterface {
        
        Q_OBJECT

    public:
        IdentificationDisplayDialog(QWidget* parent = 0);
        
        virtual ~IdentificationDisplayDialog();
        
        IdentificationDisplayDialog(const IdentificationDisplayDialog&) = delete;

        IdentificationDisplayDialog& operator=(const IdentificationDisplayDialog&) = delete;
        
        void updateDialog();
        

        // ADD_NEW_METHODS_HERE

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
        std::unique_ptr<SceneClassAssistant> m_sceneAssistant;

        IdentificationDisplayWidget* m_identificationWidget;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __IDENTIFICATION_DISPLAY_DIALOG_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __IDENTIFICATION_DISPLAY_DIALOG_DECLARE__

} // namespace
#endif  //__IDENTIFICATION_DISPLAY_DIALOG_H__
