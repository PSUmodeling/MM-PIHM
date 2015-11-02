#include "Cycles.h"

int IsOperationToday (int rotationYear, int doy, FieldOperationStruct *FieldOperation, int numOperation, int *operationIndex)
{
    /*
     * Returns a true or false indicating if an operation happens on that day
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * operation_today	    int		[return value]
     */
    int             operation_today;
    int             i, j;

    operation_today = 0;
    *operationIndex = -1;

    if  (numOperation <= 0)
    {
        operation_today = 0;
        *operationIndex = -1;
    }
    else
    {
        for (i = 0; i < numOperation; i++)
        {
            if (FieldOperation[i].status == 1)
                continue;
            else
            {
                if (rotationYear == FieldOperation[i].opYear && doy == FieldOperation[i].opDay)
                {
                    operation_today = 1;
                    FieldOperation[i].status = 1;
                    *operationIndex = i;

                    break;
                }
                else
                    break;
            }
        }
    }

    return (operation_today);
}

void UpdateOperationStatus (FieldOperationStruct *FieldOperation, int numOperation)
{
    int             i;
    int             all_performed = 1;

    if (numOperation <= 0)
    {
        /* Do nothing */
    }
    else
    {
        for (i = 0; i < numOperation; i++)
        {
            if (FieldOperation[i].status == 0)
            {
                all_performed = 0;
                break;
            }
        }

        if (all_performed)
        {
            for (i = 0 ; i < numOperation; i++)
                FieldOperation[i].status = 0;
        }
    }
}
//void SelectNextOperation (int NumOperation, int *operationIndex)
//{
//    /*
//     * Select next operation in the list, if any
//     */
//    (*operationIndex)++;
//    if (*operationIndex >= NumOperation)
//        *operationIndex = -1;
//}
//
//void SelectOperationYear (int rotationYear, const FieldOperationStruct *FieldOperation, int NumOperation, int *operationIndex)
//{
//    if (NumOperation == 0)
//    {
//        *operationIndex = -1;
//    }
//    else
//    {
//        if (rotationYear > 0)
//        {
//            *operationIndex = 0;
//            while (*operationIndex >= 0)
//            {
//                if (FieldOperation[*operationIndex].opYear < rotationYear)
//                    SelectNextOperation (NumOperation, operationIndex);
//                else
//                    break;
//            }
//        }
//        else
//        {
//            printf ("ERROR: Field operations attempting to be passed an invalid year: %d", rotationYear);
//            exit (1);
//        }
//    }
//}
//
//int IsOperationToday (int rotationYear, int doy, const FieldOperationStruct *FieldOperation, int operationIndex)
//{
//    /*
//     * Returns a true or false indicating if an operation happens on that day
//     * -----------------------------------------------------------------------
//     * LOCAL VARIABLES
//     *
//     * Variable             Type        Description
//     * ==========           ==========  ====================
//     * operation_today	    int		[return value]
//     */
//    int             operation_today;
//
//    if (operationIndex == -1)
//        operation_today = 0;
//    else
//    {
//        if (rotationYear == FieldOperation[operationIndex].opYear && doy == FieldOperation[operationIndex].opDay)
//            operation_today = 1;
//        else
//            operation_today = 0;
//    }
//
//    return (operation_today);
//}
