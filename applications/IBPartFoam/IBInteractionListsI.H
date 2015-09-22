/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::polyMesh& Foam::IBInteractionLists<ParticleType>::mesh() const
{
    return mesh_;
}


template<class ParticleType>
const Foam::mapDistribute&
Foam::IBInteractionLists<ParticleType>::cellMap() const
{
    return cellMapPtr_();
}


template<class ParticleType>
const Foam::mapDistribute&
Foam::IBInteractionLists<ParticleType>::wallFaceMap() const
{
    return wallFaceMapPtr_();
}


template<class ParticleType>
const Foam::labelListList& Foam::IBInteractionLists<ParticleType>::dil() const
{
    return dil_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::IBInteractionLists<ParticleType>::dwfil() const
{
    return dwfil_;
}


template<class ParticleType>
const Foam::labelListList& Foam::IBInteractionLists<ParticleType>::ril() const
{
    return ril_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::IBInteractionLists<ParticleType>::rilInverse() const
{
    return rilInverse_;
}


template<class ParticleType>
const Foam::labelListList& Foam::IBInteractionLists<ParticleType>::rwfil() const
{
    return rwfil_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::IBInteractionLists<ParticleType>::rwfilInverse() const
{
    return rwfilInverse_;
}


template<class ParticleType>
const Foam::List<Foam::labelPair>&
Foam::IBInteractionLists<ParticleType>::cellIndexAndTransformToDistribute() const
{
    return cellIndexAndTransformToDistribute_;
}


template<class ParticleType>
const Foam::List<Foam::labelPair>&
Foam::IBInteractionLists<ParticleType>::
wallFaceIndexAndTransformToDistribute() const
{
    return wallFaceIndexAndTransformToDistribute_;
}


template<class ParticleType>
const Foam::List<Foam::referredWallFace>&
Foam::IBInteractionLists<ParticleType>::referredWallFaces() const
{
    return referredWallFaces_;
}


template<class ParticleType>
const Foam::word& Foam::IBInteractionLists<ParticleType>::UName() const
{
    return UName_;
}


template<class ParticleType>
const Foam::List<Foam::vector>&
Foam::IBInteractionLists<ParticleType>::referredWallData() const
{
    return referredWallData_;
}


template<class ParticleType>
const Foam::List<Foam::IDLList<ParticleType> >&
Foam::IBInteractionLists<ParticleType>::referredParticles() const
{
    return referredParticles_;
}


template<class ParticleType>
Foam::List<Foam::IDLList<ParticleType> >&
Foam::IBInteractionLists<ParticleType>::referredParticles()
{
    return referredParticles_;
}


// ************************************************************************* //
