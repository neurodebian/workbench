#ifndef __WU_Q_TAB_WIDGET__H_
#define __WU_Q_TAB_WIDGET__H_

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


#include "WuQWidget.h"
#include "SceneableInterface.h"

class QTabBar;
class QStackedWidget;

namespace caret {

    class WuQTabWidget : public WuQWidget, public SceneableInterface {
        
        Q_OBJECT

    public:
        enum TabAlignment {
            TAB_ALIGN_LEFT,
            TAB_ALIGN_CENTER,
            TAB_ALIGN_RIGHT
        };
        
        WuQTabWidget(const TabAlignment alignment,
                     QObject* parent);
        
        virtual ~WuQTabWidget();
        
        QWidget* getWidget();
        
        void addTab(QWidget* page,
                    const QString& label);
        
        int currentIndex() const;
        
        QWidget* currentWidget() const;
        
        QTabBar* getTabBar() const;
        
        virtual SceneClass* saveToScene(const SceneAttributes* sceneAttributes,
                                        const AString& instanceName);
        
        virtual void restoreFromScene(const SceneAttributes* sceneAttributes,
                                      const SceneClass* sceneClass);
        
    signals:
        void currentChanged(int index);
        
    public slots:
        void setCurrentIndex(int index);
        
        void setCurrentWidget(QWidget* widget);
        
    private slots:
        void tabBarCurrentIndexChanged(int index);
        
    private:
        WuQTabWidget(const WuQTabWidget&);

        WuQTabWidget& operator=(const WuQTabWidget&);
        
        QTabBar* m_tabBar;
        
        QStackedWidget* m_stackedWidget;
        
        QWidget* m_widget;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __WU_Q_TAB_WIDGET_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __WU_Q_TAB_WIDGET_DECLARE__

} // namespace
#endif  //__WU_Q_TAB_WIDGET__H_
