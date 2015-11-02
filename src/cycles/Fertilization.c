#include "Cycles.h"

void ApplyFertilizer (FieldOperationStruct *fixedFertilization, SoilStruct *Soil, ResidueStruct *Residue)
{
    /*
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * layer		    int
     */
    int             layer;

    /* Incorporates fertilizer with option to put manure C and manure N in the
     * surface, but mineral components are added to layer 0 */
    if (fixedFertilization->opLayer > 0)
    {
        layer = fixedFertilization->opLayer - 1;
        Residue->manureC[layer] += fixedFertilization->opMass * fixedFertilization->opC_Organic;
        Residue->manureN[layer] += fixedFertilization->opMass * fixedFertilization->opN_Organic;
        Soil->NO3[layer] += fixedFertilization->opMass * fixedFertilization->opN_NO3;
        Soil->NH4[layer] += fixedFertilization->opMass * fixedFertilization->opN_NH4;
    }
    else
    {
        Residue->manureSurfaceC += fixedFertilization->opMass * fixedFertilization->opC_Organic;
        Residue->manureSurfaceN += fixedFertilization->opMass * fixedFertilization->opN_Organic;

        layer = 0;
        Soil->NO3[layer] += fixedFertilization->opMass * fixedFertilization->opN_NO3;
        Soil->NH4[layer] += fixedFertilization->opMass * fixedFertilization->opN_NH4;
    }
}
