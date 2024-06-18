
/*LICENSE_START*/
/*
 *  Copyright (C) 2015 Washington University School of Medicine
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

#include <algorithm>
#define __BRAIN_BROWSER_WINDOW_EDIT_MENU_ITEM_ENUM_DECLARE__
#include "BrainBrowserWindowEditMenuItemEnum.h"
#undef __BRAIN_BROWSER_WINDOW_EDIT_MENU_ITEM_ENUM_DECLARE__

#include "CaretAssert.h"

using namespace caret;

    
/**
 * \class caret::BrainBrowserWindowEditMenuItemEnum 
 * \brief Enumerated type for items on the Edit Menu
 *
 * Using this enumerated type in the GUI with an EnumComboBoxTemplate
 * 
 * Header File (.h)
 *     Forward declare the data type:
 *         class EnumComboBoxTemplate;
 * 
 *     Declare the member:
 *         EnumComboBoxTemplate* m_brainBrowserWindowEditMenuItemEnumComboBox;
 * 
 *     Declare a slot that is called when user changes selection
 *         private slots:
 *             void brainBrowserWindowEditMenuItemEnumComboBoxItemActivated();
 * 
 * Implementation File (.cxx)
 *     Include the header files
 *         #include "EnumComboBoxTemplate.h"
 *         #include "BrainBrowserWindowEditMenuItemEnum.h"
 * 
 *     Instatiate:
 *         m_brainBrowserWindowEditMenuItemEnumComboBox = new EnumComboBoxTemplate(this);
 *         m_brainBrowserWindowEditMenuItemEnumComboBox->setup<BrainBrowserWindowEditMenuItemEnum,BrainBrowserWindowEditMenuItemEnum::Enum>();
 * 
 *     Get notified when the user changes the selection: 
 *         QObject::connect(m_brainBrowserWindowEditMenuItemEnumComboBox, SIGNAL(itemActivated()),
 *                          this, SLOT(brainBrowserWindowEditMenuItemEnumComboBoxItemActivated()));
 * 
 *     Update the selection:
 *         m_brainBrowserWindowEditMenuItemEnumComboBox->setSelectedItem<BrainBrowserWindowEditMenuItemEnum,BrainBrowserWindowEditMenuItemEnum::Enum>(NEW_VALUE);
 * 
 *     Read the selection:
 *         const BrainBrowserWindowEditMenuItemEnum::Enum VARIABLE = m_brainBrowserWindowEditMenuItemEnumComboBox->getSelectedItem<BrainBrowserWindowEditMenuItemEnum,BrainBrowserWindowEditMenuItemEnum::Enum>();
 * 
 */

/**
 * Constructor.
 *
 * @param enumValue
 *    An enumerated value.
 * @param name
 *    Name of enumerated value.
 * @param guiName
 *    User-friendly name for use in user-interface.
 * @param shortCut
 *    Shortcut for menu item.
 */
BrainBrowserWindowEditMenuItemEnum::BrainBrowserWindowEditMenuItemEnum(const Enum enumValue,
                                                                       const AString& name,
                                                                       const AString& guiName,
                                                                       const QKeySequence& shortCut)
{
    this->enumValue   = enumValue;
    this->integerCode = integerCodeCounter++;
    this->name        = name;
    this->guiName     = guiName;
    this->shortCut    = shortCut;
}

/**
 * Destructor.
 */
BrainBrowserWindowEditMenuItemEnum::~BrainBrowserWindowEditMenuItemEnum()
{
}

/**
 * Initialize the enumerated metadata.
 */
void
BrainBrowserWindowEditMenuItemEnum::initialize()
{
    if (initializedFlag) {
        return;
    }
    initializedFlag = true;

    QKeySequence noShortCutKeySequence;

    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(COPY,
                                                          "COPY",
                                                          "Copy",
                                                          QKeySequence(Qt::CTRL | Qt::Key_C)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(CUT,
                                                          "CUT",
                                                          "Cut",
                                                          QKeySequence(Qt::CTRL | Qt::Key_X)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(DELETER,
                                                          "DELETER",
                                                          "Delete",
                                                          noShortCutKeySequence));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(DESELECT_ALL,
                                                          "DESELECT_ALL",
                                                          "Deselect All",
                                                          QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_A)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(PASTE,
                                                          "PASTE",
                                                          "Paste",
                                                          QKeySequence(Qt::CTRL | Qt::Key_V)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(PASTE_SPECIAL,
                                                          "PASTE",
                                                          "Paste Special",
                                                          QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_V)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(REDO,
                                                          "REDO",
                                                          "Redo",
                                                          QKeySequence(Qt::CTRL | Qt::Key_Y)));  //(Qt::CTRL + Qt::SHIFT + Qt::Key_Z)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(SELECT_ALL,
                                                          "SELECT_ALL",
                                                          "Select All",
                                                          QKeySequence(Qt::CTRL | Qt::Key_A)));
    
    enumData.push_back(BrainBrowserWindowEditMenuItemEnum(UNDO, 
                                                          "UNDO", 
                                                          "Undo",
                                                          QKeySequence(Qt::CTRL | Qt::Key_Z)));
}

/**
 * Find the data for and enumerated value.
 * @param enumValue
 *     The enumerated value.
 * @return Pointer to data for this enumerated type
 * or NULL if no data for type or if type is invalid.
 */
const BrainBrowserWindowEditMenuItemEnum*
BrainBrowserWindowEditMenuItemEnum::findData(const Enum enumValue)
{
    if (initializedFlag == false) initialize();

    size_t num = enumData.size();
    for (size_t i = 0; i < num; i++) {
        const BrainBrowserWindowEditMenuItemEnum* d = &enumData[i];
        if (d->enumValue == enumValue) {
            return d;
        }
    }

    return NULL;
}

/**
 * Get a string representation of the enumerated type.
 * @param enumValue 
 *     Enumerated value.
 * @return 
 *     String representing enumerated value.
 */
AString 
BrainBrowserWindowEditMenuItemEnum::toName(Enum enumValue) {
    if (initializedFlag == false) initialize();
    
    const BrainBrowserWindowEditMenuItemEnum* enumInstance = findData(enumValue);
    return enumInstance->name;
}

/**
 * Get an enumerated value corresponding to its name.
 * @param name 
 *     Name of enumerated value.
 * @param isValidOut 
 *     If not NULL, it is set indicating that a
 *     enum value exists for the input name.
 * @return 
 *     Enumerated value.
 */
BrainBrowserWindowEditMenuItemEnum::Enum 
BrainBrowserWindowEditMenuItemEnum::fromName(const AString& name, bool* isValidOut)
{
    if (initializedFlag == false) initialize();
    
    bool validFlag = false;
    Enum enumValue = BrainBrowserWindowEditMenuItemEnum::enumData[0].enumValue;
    
    for (std::vector<BrainBrowserWindowEditMenuItemEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        const BrainBrowserWindowEditMenuItemEnum& d = *iter;
        if (d.name == name) {
            enumValue = d.enumValue;
            validFlag = true;
            break;
        }
    }
    
    if (isValidOut != 0) {
        *isValidOut = validFlag;
    }
    else if (validFlag == false) {
        CaretAssertMessage(0, AString("Name " + name + "failed to match enumerated value for type BrainBrowserWindowEditMenuItemEnum"));
    }
    return enumValue;
}

/**
 * Get a GUI string representation of the enumerated type.
 * @param enumValue 
 *     Enumerated value.
 * @return 
 *     String representing enumerated value.
 */
AString 
BrainBrowserWindowEditMenuItemEnum::toGuiName(Enum enumValue) {
    if (initializedFlag == false) initialize();
    
    const BrainBrowserWindowEditMenuItemEnum* enumInstance = findData(enumValue);
    return enumInstance->guiName;
}

/**
 * Get an enumerated value corresponding to its GUI name.
 * @param s 
 *     Name of enumerated value.
 * @param isValidOut 
 *     If not NULL, it is set indicating that a
 *     enum value exists for the input name.
 * @return 
 *     Enumerated value.
 */
BrainBrowserWindowEditMenuItemEnum::Enum 
BrainBrowserWindowEditMenuItemEnum::fromGuiName(const AString& guiName, bool* isValidOut)
{
    if (initializedFlag == false) initialize();
    
    bool validFlag = false;
    Enum enumValue = BrainBrowserWindowEditMenuItemEnum::enumData[0].enumValue;
    
    for (std::vector<BrainBrowserWindowEditMenuItemEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        const BrainBrowserWindowEditMenuItemEnum& d = *iter;
        if (d.guiName == guiName) {
            enumValue = d.enumValue;
            validFlag = true;
            break;
        }
    }
    
    if (isValidOut != 0) {
        *isValidOut = validFlag;
    }
    else if (validFlag == false) {
        CaretAssertMessage(0, AString("guiName " + guiName + "failed to match enumerated value for type BrainBrowserWindowEditMenuItemEnum"));
    }
    return enumValue;
}

/**
 * Get the integer code for a data type.
 *
 * @return
 *    Integer code for data type.
 */
int32_t
BrainBrowserWindowEditMenuItemEnum::toIntegerCode(Enum enumValue)
{
    if (initializedFlag == false) initialize();
    const BrainBrowserWindowEditMenuItemEnum* enumInstance = findData(enumValue);
    return enumInstance->integerCode;
}

/**
 * Find the data type corresponding to an integer code.
 *
 * @param integerCode
 *     Integer code for enum.
 * @param isValidOut
 *     If not NULL, on exit isValidOut will indicate if
 *     integer code is valid.
 * @return
 *     Enum for integer code.
 */
BrainBrowserWindowEditMenuItemEnum::Enum
BrainBrowserWindowEditMenuItemEnum::fromIntegerCode(const int32_t integerCode, bool* isValidOut)
{
    if (initializedFlag == false) initialize();
    
    bool validFlag = false;
    Enum enumValue = BrainBrowserWindowEditMenuItemEnum::enumData[0].enumValue;
    
    for (std::vector<BrainBrowserWindowEditMenuItemEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        const BrainBrowserWindowEditMenuItemEnum& enumInstance = *iter;
        if (enumInstance.integerCode == integerCode) {
            enumValue = enumInstance.enumValue;
            validFlag = true;
            break;
        }
    }
    
    if (isValidOut != 0) {
        *isValidOut = validFlag;
    }
    else if (validFlag == false) {
        CaretAssertMessage(0, AString("Integer code " + AString::number(integerCode) + "failed to match enumerated value for type BrainBrowserWindowEditMenuItemEnum"));
    }
    return enumValue;
}

/**
 * Get the shortcut key for a data type.
 *
 * @return
 *    Shortcut keys for data type.
 */
QKeySequence
BrainBrowserWindowEditMenuItemEnum::toShortcut(Enum enumValue)
{
    if (initializedFlag == false) initialize();
    const BrainBrowserWindowEditMenuItemEnum* enumInstance = findData(enumValue);
    return enumInstance->shortCut;
}

/**
 * Get all of the enumerated type values.  The values can be used
 * as parameters to toXXX() methods to get associated metadata.
 *
 * @param allEnums
 *     A vector that is OUTPUT containing all of the enumerated values.
 */
void
BrainBrowserWindowEditMenuItemEnum::getAllEnums(std::vector<BrainBrowserWindowEditMenuItemEnum::Enum>& allEnums)
{
    if (initializedFlag == false) initialize();
    
    allEnums.clear();
    
    for (std::vector<BrainBrowserWindowEditMenuItemEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        allEnums.push_back(iter->enumValue);
    }
}

/**
 * Get all of the names of the enumerated type values.
 *
 * @param allNames
 *     A vector that is OUTPUT containing all of the names of the enumerated values.
 * @param isSorted
 *     If true, the names are sorted in alphabetical order.
 */
void
BrainBrowserWindowEditMenuItemEnum::getAllNames(std::vector<AString>& allNames, const bool isSorted)
{
    if (initializedFlag == false) initialize();
    
    allNames.clear();
    
    for (std::vector<BrainBrowserWindowEditMenuItemEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        allNames.push_back(BrainBrowserWindowEditMenuItemEnum::toName(iter->enumValue));
    }
    
    if (isSorted) {
        std::sort(allNames.begin(), allNames.end());
    }
}

/**
 * Get all of the GUI names of the enumerated type values.
 *
 * @param allNames
 *     A vector that is OUTPUT containing all of the GUI names of the enumerated values.
 * @param isSorted
 *     If true, the names are sorted in alphabetical order.
 */
void
BrainBrowserWindowEditMenuItemEnum::getAllGuiNames(std::vector<AString>& allGuiNames, const bool isSorted)
{
    if (initializedFlag == false) initialize();
    
    allGuiNames.clear();
    
    for (std::vector<BrainBrowserWindowEditMenuItemEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        allGuiNames.push_back(BrainBrowserWindowEditMenuItemEnum::toGuiName(iter->enumValue));
    }
    
    if (isSorted) {
        std::sort(allGuiNames.begin(), allGuiNames.end());
    }
}

