
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

#include <algorithm>
#define __CZI_PIXEL_COORD_SPACE_ENUM_DECLARE__
#include "CziPixelCoordSpaceEnum.h"
#undef __CZI_PIXEL_COORD_SPACE_ENUM_DECLARE__

#include "CaretAssert.h"

using namespace caret;

    
/**
 * \class caret::CziPixelCoordSpaceEnum
 * \brief Enumerated type for a pixel coordinate type in CZI Image
 *
 * Using this enumerated type in the GUI with an EnumComboBoxTemplate
 * 
 * Header File (.h)
 *     Forward declare the data type:
 *         class EnumComboBoxTemplate;
 * 
 *     Declare the member:
 *         EnumComboBoxTemplate* m_cziPixelCoordTypeEnumComboBox;
 * 
 *     Declare a slot that is called when user changes selection
 *         private slots:
 *             void cziPixelCoordTypeEnumComboBoxItemActivated();
 * 
 * Implementation File (.cxx)
 *     Include the header files
 *         #include "EnumComboBoxTemplate.h"
 *         #include "CziPixelCoordSpaceEnum.h"
 * 
 *     Instatiate:
 *         m_cziPixelCoordTypeEnumComboBox = new EnumComboBoxTemplate(this);
 *         m_cziPixelCoordTypeEnumComboBox->setup<CziPixelCoordSpaceEnum,CziPixelCoordSpaceEnum::Enum>();
 * 
 *     Get notified when the user changes the selection: 
 *         QObject::connect(m_cziPixelCoordTypeEnumComboBox, SIGNAL(itemActivated()),
 *                          this, SLOT(cziPixelCoordTypeEnumComboBoxItemActivated()));
 * 
 *     Update the selection:
 *         m_cziPixelCoordTypeEnumComboBox->setSelectedItem<CziPixelCoordSpaceEnum,CziPixelCoordSpaceEnum::Enum>(NEW_VALUE);
 * 
 *     Read the selection:
 *         const CziPixelCoordSpaceEnum::Enum VARIABLE = m_cziPixelCoordTypeEnumComboBox->getSelectedItem<CziPixelCoordSpaceEnum,CziPixelCoordSpaceEnum::Enum>();
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
 * @param description
 *    Description of this coordinate type
 */
CziPixelCoordSpaceEnum::CziPixelCoordSpaceEnum(const Enum enumValue,
                                             const AString& name,
                                             const AString& guiName,
                                             const AString& description)
{
    this->enumValue = enumValue;
    this->integerCode = integerCodeCounter++;
    this->name = name;
    this->guiName = guiName;
    this->description = description;
}

/**
 * Destructor.
 */
CziPixelCoordSpaceEnum::~CziPixelCoordSpaceEnum()
{
}

/**
 * Initialize the enumerated metadata.
 */
void
CziPixelCoordSpaceEnum::initialize()
{
    if (initializedFlag) {
        return;
    }
    enumData.push_back(CziPixelCoordSpaceEnum(FULL_RESOLUTION_LOGICAL_TOP_LEFT,
                                              "FULL_RESOLUTION_LOGICAL_TOP_LEFT",
                                              "FULL_RESOLUTION_LOGICAL_TOP_LEFT Top Left",
                                              "Logical Coordinate of full resolution image with origin (x, y) at top left"));
    
    enumData.push_back(CziPixelCoordSpaceEnum(FULL_RESOLUTION_PIXEL_BOTTOM_LEFT,
                                              "FULL_RESOLUTION_PIXEL_BOTTOM_LEFT",
                                              "Full Resolution Pixel Bottom Left",
                                              "Pixel coordinate of full-resolution image with origin (0, 0) at bottom left"));
    
    enumData.push_back(CziPixelCoordSpaceEnum(FULL_RESOLUTION_PIXEL_TOP_LEFT,
                                              "FULL_RESOLUTION_PIXEL_TOP_LEFT",
                                              "Full Resolution Pixel Top Left",
                                              "Pixel coordinate of full-resolution image with origin (0, 0) at top left"));
    
    enumData.push_back(CziPixelCoordSpaceEnum(IMAGE_DATA_PIXEL_BOTTOM_LEFT,
                                             "IMAGE_DATA_PIXEL_BOTTOM_LEFT",
                                             "Image Data Pixel Bottom Left",
                                             "Pixel coordinate of image with origin (0,0) at the bottom left"));
    
    enumData.push_back(CziPixelCoordSpaceEnum(IMAGE_DATA_PIXEL_TOP_LEFT,
                                             "IMAGE_DATA_PIXEL_TOP_LEFT",
                                             "Image Data Pixel Top Left",
                                             "Pixel coordinate of image with origin (0,0) at the top left"));
}

/**
 * Find the data for and enumerated value.
 * @param enumValue
 *     The enumerated value.
 * @return Pointer to data for this enumerated type
 * or NULL if no data for type or if type is invalid.
 */
const CziPixelCoordSpaceEnum*
CziPixelCoordSpaceEnum::findData(const Enum enumValue)
{
    if (initializedFlag == false) initialize();

    size_t num = enumData.size();
    for (size_t i = 0; i < num; i++) {
        const CziPixelCoordSpaceEnum* d = &enumData[i];
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
CziPixelCoordSpaceEnum::toName(Enum enumValue) {
    if (initializedFlag == false) initialize();
    
    const CziPixelCoordSpaceEnum* enumInstance = findData(enumValue);
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
CziPixelCoordSpaceEnum::Enum
CziPixelCoordSpaceEnum::fromName(const AString& name, bool* isValidOut)
{
    if (initializedFlag == false) initialize();
    
    bool validFlag = false;
    Enum enumValue = CziPixelCoordSpaceEnum::enumData[0].enumValue;
    
    for (std::vector<CziPixelCoordSpaceEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        const CziPixelCoordSpaceEnum& d = *iter;
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
        CaretAssertMessage(0, AString("Name " + name + " failed to match enumerated value for type CziPixelCoordSpaceEnum"));
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
CziPixelCoordSpaceEnum::toGuiName(Enum enumValue) {
    if (initializedFlag == false) initialize();
    
    const CziPixelCoordSpaceEnum* enumInstance = findData(enumValue);
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
CziPixelCoordSpaceEnum::Enum
CziPixelCoordSpaceEnum::fromGuiName(const AString& guiName, bool* isValidOut)
{
    if (initializedFlag == false) initialize();
    
    bool validFlag = false;
    Enum enumValue = CziPixelCoordSpaceEnum::enumData[0].enumValue;
    
    for (std::vector<CziPixelCoordSpaceEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        const CziPixelCoordSpaceEnum& d = *iter;
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
        CaretAssertMessage(0, AString("guiName " + guiName + " failed to match enumerated value for type CziPixelCoordSpaceEnum"));
    }
    return enumValue;
}

/**
 * Get a GUI string representation of the description.
 * @param enumValue
 *     Enumerated value.
 * @return
 *     String representing description.
 */
AString
CziPixelCoordSpaceEnum::toDescription(Enum enumValue) {
    if (initializedFlag == false) initialize();
    
    const CziPixelCoordSpaceEnum* enumInstance = findData(enumValue);
    return enumInstance->description;
}

/**
 * Get the integer code for a data type.
 *
 * @return
 *    Integer code for data type.
 */
int32_t
CziPixelCoordSpaceEnum::toIntegerCode(Enum enumValue)
{
    if (initializedFlag == false) initialize();
    const CziPixelCoordSpaceEnum* enumInstance = findData(enumValue);
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
CziPixelCoordSpaceEnum::Enum
CziPixelCoordSpaceEnum::fromIntegerCode(const int32_t integerCode, bool* isValidOut)
{
    if (initializedFlag == false) initialize();
    
    bool validFlag = false;
    Enum enumValue = CziPixelCoordSpaceEnum::enumData[0].enumValue;
    
    for (std::vector<CziPixelCoordSpaceEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        const CziPixelCoordSpaceEnum& enumInstance = *iter;
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
        CaretAssertMessage(0, AString("Integer code " + AString::number(integerCode) + " failed to match enumerated value for type CziPixelCoordSpaceEnum"));
    }
    return enumValue;
}

/**
 * Get all of the enumerated type values.  The values can be used
 * as parameters to toXXX() methods to get associated metadata.
 *
 * @param allEnums
 *     A vector that is OUTPUT containing all of the enumerated values.
 */
void
CziPixelCoordSpaceEnum::getAllEnums(std::vector<CziPixelCoordSpaceEnum::Enum>& allEnums)
{
    if (initializedFlag == false) initialize();
    
    allEnums.clear();
    
    for (std::vector<CziPixelCoordSpaceEnum>::iterator iter = enumData.begin();
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
CziPixelCoordSpaceEnum::getAllNames(std::vector<AString>& allNames, const bool isSorted)
{
    if (initializedFlag == false) initialize();
    
    allNames.clear();
    
    for (std::vector<CziPixelCoordSpaceEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        allNames.push_back(CziPixelCoordSpaceEnum::toName(iter->enumValue));
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
CziPixelCoordSpaceEnum::getAllGuiNames(std::vector<AString>& allGuiNames, const bool isSorted)
{
    if (initializedFlag == false) initialize();
    
    allGuiNames.clear();
    
    for (std::vector<CziPixelCoordSpaceEnum>::iterator iter = enumData.begin();
         iter != enumData.end();
         iter++) {
        allGuiNames.push_back(CziPixelCoordSpaceEnum::toGuiName(iter->enumValue));
    }
    
    if (isSorted) {
        std::sort(allGuiNames.begin(), allGuiNames.end());
    }
}

