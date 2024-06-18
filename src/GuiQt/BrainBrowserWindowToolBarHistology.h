#ifndef __BRAIN_BROWSER_WINDOW_TOOL_BAR_HISTOLOGY_H__
#define __BRAIN_BROWSER_WINDOW_TOOL_BAR_HISTOLOGY_H__

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

#include "BrainBrowserWindowToolBarComponent.h"

class QCheckBox;
class QComboBox;
class QLabel;

namespace caret {

    class BrainOpenGLViewportContent;
    class BrowserTabContent;
    class HistologySlice;
    class HistologySlicesFile;
    class Vector3D;
    class WuQDoubleSpinBox;
    class WuQSpinBox;

    class BrainBrowserWindowToolBarHistology : public BrainBrowserWindowToolBarComponent {
        Q_OBJECT
        
    public:
        BrainBrowserWindowToolBarHistology(BrainBrowserWindowToolBar* parentToolBar,
                                           const QString& parentObjectName);
        
        virtual ~BrainBrowserWindowToolBarHistology();
        
        virtual void updateContent(BrowserTabContent* browserTabContent);
        
        void receiveEvent(Event* event) override;
        
    private:
        BrainBrowserWindowToolBarHistology(const BrainBrowserWindowToolBarHistology&);

        BrainBrowserWindowToolBarHistology& operator=(const BrainBrowserWindowToolBarHistology&);
        
    public:

        // ADD_NEW_METHODS_HERE

    private slots:
        void sliceIndexValueChanged(int);
        
        void sliceNameComboBoxActivated(int);
        
        void planeXyzSpinBoxValueChanged();

        void stereotaxicXyzSpinBoxValueChanged();
        
        void identificationMovesSlicesActionTriggered(bool);
        
        void moveToCenterActionTriggered();
        
        void yokeOrientationCheckBoxChecked(bool checked);
        
        void axisCrosshairActionTriggered(bool checked);
        
    private:
        HistologySlicesFile* getHistologySlicesFile(BrowserTabContent* browserTabContent);
        
        const BrainOpenGLViewportContent* getBrainOpenGLViewportContent() const;
        
        bool getPlaneCoordinateAtViewportCenter(Vector3D& planeXyzOut) const;
        
        bool getStereotaxicCoordinateAtViewportCenter(const HistologySlice* histologySlice,
                                                      Vector3D& stereotaxicXyzOut) const;
        // ADD_NEW_MEMBERS_HERE

        BrainBrowserWindowToolBar* m_parentToolBar;
        
        BrowserTabContent* m_browserTabContent = NULL;
        
        WuQSpinBox* m_sliceIndexSpinBox;
        
        QComboBox* m_sliceNameComboBox;
        
        WuQDoubleSpinBox* m_planeXyzSpinBox[2];
        
        WuQDoubleSpinBox* m_stereotaxicXyzSpinBox[3];
        
        QAction* m_identificationMovesSlicesAction;
        
        QAction* m_moveToCenterAction;
        
        QCheckBox* m_yokeOrientationCheckBox;
        
        QAction* m_showAxisCrosshairsAction;
        
        QLabel* m_rotationAngleXLabel;
        
        QLabel* m_rotationAngleYLabel;
        
        QLabel* m_rotationAngleZLabel;
        
};
    
    
#ifdef __BRAIN_BROWSER_WINDOW_TOOL_BAR_HISTOLOGY_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __BRAIN_BROWSER_WINDOW_TOOL_BAR_HISTOLOGY_DECLARE__

} // namespace
#endif  //__BRAIN_BROWSER_WINDOW_TOOL_BAR_HISTOLOGY_H__
