#ifndef __WU_Q_VALUE_CHANGED_SIGNAL_WATCHER_H__
#define __WU_Q_VALUE_CHANGED_SIGNAL_WATCHER_H__

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
#include <vector>

#include <QObject>

class QWidget;

namespace caret {

    class WuQValueChangedSignalWatcher : public QObject {
        
        Q_OBJECT

    public:
        WuQValueChangedSignalWatcher(QObject* parent);
        
        virtual ~WuQValueChangedSignalWatcher();
        
        WuQValueChangedSignalWatcher(const WuQValueChangedSignalWatcher&) = delete;

        WuQValueChangedSignalWatcher& operator=(const WuQValueChangedSignalWatcher&) = delete;
        
        void addObject(QObject* object);
        
        void setWidgetsVisible(const bool status);
        
    signals:
        void valueChanged();

//        private slots:
//        void booChanged(bool);
//        
//        void intChanged(int);
//        
//        void doubleChanged(double);
        
        // ADD_NEW_METHODS_HERE

    private:
        std::vector<QWidget*> m_widgets;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __WU_Q_VALUE_CHANGED_SIGNAL_WATCHER_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __WU_Q_VALUE_CHANGED_SIGNAL_WATCHER_DECLARE__

} // namespace
#endif  //__WU_Q_VALUE_CHANGED_SIGNAL_WATCHER_H__
