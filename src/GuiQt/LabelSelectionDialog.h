#ifndef __LABEL_SELECTION_DIALOG_H__
#define __LABEL_SELECTION_DIALOG_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2023 Washington University School of Medicine
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

#include "WuQDialogModal.h"

namespace caret {

    class CaretMappableDataFile;
    class LabelSelectionWidget;
    
    class LabelSelectionDialog : public WuQDialogModal {
        
        Q_OBJECT

    public:
        enum class Mode {
            FILE_AND_MAP,
            FILE_MAP_AND_LABEL
        };

        LabelSelectionDialog(const Mode mode,
                             const QString& saveRestoreStateName = "",
                             QWidget* parent = 0);
        
        virtual ~LabelSelectionDialog();
        
        LabelSelectionDialog(const LabelSelectionDialog&) = delete;

        LabelSelectionDialog& operator=(const LabelSelectionDialog&) = delete;
        
        AString getSelectedLabel() const;

        CaretMappableDataFile* getSelectedFile() const;
        
        AString getSelectedFileName() const;
        
        AString getSelectedFileNameNoPath() const;
        
        AString getSelectedMapName() const;
        
        int32_t getSelectedMapIndex() const;

        // ADD_NEW_METHODS_HERE

    private:
        LabelSelectionWidget* m_labelSelectionWidget = NULL;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __LABEL_SELECTION_DIALOG_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __LABEL_SELECTION_DIALOG_DECLARE__

} // namespace
#endif  //__LABEL_SELECTION_DIALOG_H__
