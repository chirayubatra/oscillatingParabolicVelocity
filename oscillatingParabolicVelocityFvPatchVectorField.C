/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 Chirayu Batra
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "oscillatingParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oscillatingParabolicVelocityFvPatchVectorField::
oscillatingParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0.0),
    n_(1, 0, 0),
    y_(0, 1, 0), 
    freq_(0.0)
{
}


Foam::oscillatingParabolicVelocityFvPatchVectorField::
oscillatingParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(dict.get<scalar>("maxValue")),
    n_(dict.get<vector>("n")),
    y_(dict.get<vector>("y")),
    freq_(dict.get<scalar>("freq"))
{
    Info << "Using the oscillatingParabolicVelocity boundary condition" << endl;
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("oscillatingParabolicVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }
    n_ /= mag(n_);
    y_ /= mag(y_);

    fixedValueFvPatchVectorField::evaluate();

    /*
    //Initialise with the value entry if evaluation is not possible
    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
    */
}


Foam::oscillatingParabolicVelocityFvPatchVectorField::
oscillatingParabolicVelocityFvPatchVectorField
(
    const oscillatingParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_), 
    freq_(ptf.freq_)
{}


Foam::oscillatingParabolicVelocityFvPatchVectorField::
oscillatingParabolicVelocityFvPatchVectorField
(
    const oscillatingParabolicVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_), 
    freq_(ptf.freq_)
{}


Foam::oscillatingParabolicVelocityFvPatchVectorField::
oscillatingParabolicVelocityFvPatchVectorField
(
    const oscillatingParabolicVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_),
    freq_(ptf.freq_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::oscillatingParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Get range and orientation
    boundBox bb(patch().patch().localPoints(), true);

    vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    const scalar pi = Foam::constant::mathematical::pi;
    const scalar t = this->db().time().timeOutputValue(); 
//    const scalar t = db().time().timeOutputValue(); 
//    const scalar t = mesh_.time().timeOutputValue(); 
//    const scalar frequency = freqPtr_->value(t); 

    // Calculate local 1-D coordinate for the parabolic profile
    scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);
    scalar oscMult = fabs(sin(2.0*pi*freq_*t));  //get the absolute(fabs) value 
//    coord = coord * oscMult; 
    vectorField::operator=(n_*maxValue_*(1.0 - sqr(coord))*oscMult);

}


void Foam::oscillatingParabolicVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
//    fvPatchVectorField::write(os);
//    os.writeEntry("scalarData", scalarData_);
//    os.writeEntry("data", data_);
//    fieldData_.writeEntry("fieldData", os);
//    timeVsData_->writeData(os);
//    os.writeEntry("wordData", wordData_);
//    writeEntry("value", os);
    
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue") << maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y") << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("freq") << freq_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        oscillatingParabolicVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
