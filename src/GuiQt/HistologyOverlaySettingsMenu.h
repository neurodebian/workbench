#ifndef __HISTOLOGY_OVERLAY_SETTINGS_MENU_H__
#define __HISTOLOGY_OVERLAY_SETTINGS_MENU_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2021 Washington University School of Medicine
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

#include <QMenu>

class QSpinBox;


namespace caret {

    class EnumComboBoxTemplate;
    class HistologyOverlay;
    
    class HistologyOverlaySettingsMenu : public QMenu {
        
        Q_OBJECT

    public:
        HistologyOverlaySettingsMenu(HistologyOverlay* histologyOverlay,
                                     const QString& parentObjectName);
        
        virtual ~HistologyOverlaySettingsMenu();
        
        HistologyOverlaySettingsMenu(const HistologyOverlaySettingsMenu&) = delete;

        HistologyOverlaySettingsMenu& operator=(const HistologyOverlaySettingsMenu&) = delete;
        

        // ADD_NEW_METHODS_HERE

    private slots:
        void reloadActionTriggered();
        
    private:
        void updateContent();
        
        HistologyOverlay* m_histologyOverlay;
        
        QAction* m_reloadAction;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __HISTOLOGY_OVERLAY_SETTINGS_MENU_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __HISTOLOGY_OVERLAY_SETTINGS_MENU_DECLARE__

} // namespace
#endif  //__HISTOLOGY_OVERLAY_SETTINGS_MENU_H__
