#ifndef __META_DATA_CUSTOM_EDITOR_DIALOG_H__
#define __META_DATA_CUSTOM_EDITOR_DIALOG_H__

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


#include "WuQDialogModal.h"

namespace caret {

    class Annotation;
    class CaretDataFile;
    class CaretMappableDataFile;
    class GiftiMetaData;
    class MetaDataCustomEditorWidget;
    
    class MetaDataCustomEditorDialog : public WuQDialogModal {
        
        Q_OBJECT

    public:
        MetaDataCustomEditorDialog(Annotation* annotation,
                             QWidget* parent);
        
        MetaDataCustomEditorDialog(CaretDataFile* caretDataFile,
                             QWidget* parent);
        
        MetaDataCustomEditorDialog(CaretMappableDataFile* caretMappableDataFile,
                             const int32_t mapIndex,
                             QWidget* parent);
        
        virtual ~MetaDataCustomEditorDialog();
        
        virtual void okButtonClicked();
        
        virtual void cancelButtonClicked();
        
    private:
        MetaDataCustomEditorDialog(const MetaDataCustomEditorDialog&);

        MetaDataCustomEditorDialog& operator=(const MetaDataCustomEditorDialog&);
        
        void initializeDialog(const AString& dialogTitle,
                              const std::vector<AString>& metaDataNames,
                              GiftiMetaData* metaData);
        
        MetaDataCustomEditorWidget* m_metaDataEditorWidget;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __META_DATA_CUSTOM_EDITOR_DIALOG_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __META_DATA_CUSTOM_EDITOR_DIALOG_DECLARE__

} // namespace
#endif  //__META_DATA_CUSTOM_EDITOR_DIALOG_H__
