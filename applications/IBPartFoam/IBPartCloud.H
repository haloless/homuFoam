/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::IBPartCloud

Description

SourceFiles
    IBPartCloudI.H
    IBPartCloud.C
    IBPartCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef IBPartCloud_H
#define IBPartCloud_H

#include "KinematicCloud.H"
#include "IBInteractionLists.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



/*---------------------------------------------------------------------------*\
                         Class IBPartCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IBPartCloud
:
    public CloudType
{
    

public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef IBPartCloud<CloudType> ibPartCloudType;

private:
    // Private data

        //- Cloud copy pointer
        autoPtr<IBPartCloud<CloudType> > cloudCopyPtr_;
        
        //- Interactions lists determining which cells are in
        //  interaction range of each other
        IBInteractionLists<typename CloudType::parcelType> ibil_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IBPartCloud(const IBPartCloud&);

        //- Disallow default bitwise assignment
        void operator=(const IBPartCloud&);

protected:
    // Protected data
    
        // IB properties
        //dictionary ibModelProperties_;
        
        
        // Cloud evolution functions
            //- 
            
            //- Reset state of cloud
            void cloudReset(IBPartCloud<CloudType>& c)
            {
                CloudType::cloudReset(c);
            }
            
            //- 
            scalar calcPartInCellFrac(label cellJ, 
                const vector &xp, scalar rp) const
            {
                const scalar dh = pow(this->mesh().V()[cellJ], 1.0/3.0);
                const scalar thickness = dh * 0.5;
                const vector xc = this->mesh().C()[cellJ];
                
                //rp *= 10; // TODO
                scalar dist = mag(xc - xp) - rp;
                
                scalar alpha = 0;
                if (dist > thickness) {
                    alpha = 0;
                } else if (dist < -thickness) {
                    alpha = 1;
                } else {
                    alpha = 0.5 - dist/thickness*0.5;
                }
                
                return alpha;
            }
            
    
    
public:
    // Static data members

        //- Static data staticData
        //static const dataType staticData;


    // Constructors

        //- Construct given carrier gas fields
        IBPartCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        )
        :
            CloudType(cloudName, rho, U, mu, g, false),
            //ibModelProperties_(this->particleProperties().subDict("ibModel")),
            ibil_(
                this->mesh_, 
                readScalar(this->ibModelProperties().lookup("maxIBRadius")), 
                false,
                U.name()
            )
        {
            if (this->solution().steadyState())
            {
                FatalErrorIn("IBPartCloud") << "steady state not supported" 
                    << exit(FatalError);
            }
            
            if (this->solution().active()) 
            {
                if (readFields)
                {
                    parcelType::readFields(*this);
                }
            }
        }

        //- Copy constructor with new name
        IBPartCloud
        (
            IBPartCloud<CloudType>& c,
            const word& name
        )
        :
            CloudType(c, name),
            //ibModelProperties_(c.ibModelProperties_),
            ibil_(
                c.mesh(), 
                readScalar(this->ibModelProperties().lookup("maxIBRadius")), 
                false,
                //U.name()
                "U.air" // TODO
            )
        {}

        //- Copy constructor with new name - creates bare cloud
        IBPartCloud
        (
            const fvMesh& mesh,
            const word& name,
            const IBPartCloud<CloudType>& c
        )
        :
            CloudType(mesh, name, c),
            //ibModelProperties_(c.ibModelProperties_),
            ibil_(
                mesh,
                readScalar(this->ibModelProperties().lookup("maxIBRadius")), 
                false,
                //U.name()
                "U.air" // TODO
            )
        {}
        
        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new IBPartCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new IBPartCloud(this->mesh(), name, *this)
            );
        }



    // Selectors

        //- Select null constructed
        //static autoPtr<IBPartCloud> New();


    //- Destructor
    virtual ~IBPartCloud()
    {}


    // Member Functions
    
        // Access
        
            //- Return a reference to the cloud copy
            inline const IBPartCloud& cloudCopy() const {
                return cloudCopyPtr_();
            }
            
            //- 
            inline const dictionary& ibModelProperties() const {
                //return ibModelProperties_;
                return this->particleProperties().subDict("ibModel");
            }

        // Check
        
            // Fields
                //
                inline tmp<volScalarField> ibfrac() 
                {
                    tmp<volScalarField> tmp_ibfrac
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                this->name() + ":ibfrac",
                                this->db().time().timeName(),
                                this->db(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE,
                                false
                            ),
                            this->mesh_,
                            dimensionedScalar("zero", dimless, 0.0),
                            zeroGradientFvPatchScalarField::typeName
                        )
                    );
                    
                    volScalarField &frac = tmp_ibfrac();
                    
                    //
                    PstreamBuffers pBufs(Pstream::nonBlocking);
                    label startOfRequests = Pstream::nRequests();
                    ibil_.sendReferredData(this->cellOccupancy(), pBufs);
                    
                    // real-real
                    { 
                        // direct interaction list (DIL) 
                        const labelListList& dil = ibil_.dil();
                        
                        List<DynamicList<parcelType*> > &cellOccupancy = 
                            this->cellOccupancy();
                        
                        forAll(dil, realCellI)
                        {
                            // loop all particles in Cell A
                            forAll(cellOccupancy[realCellI], a)
                            {
                                parcelType *pA = cellOccupancy[realCellI][a];
                                const vector xp = pA->position();
                                const scalar dp = pA->d();
                                const scalar rp = dp * 0.5;
                                
                                // Cell B interacting with A
                                forAll(dil[realCellI], interactingCells)
                                {
                                    //List<parcelType*> cellBParcels = 
                                    //    cellOccupancy[dil[realCellI][interactingCells]];
                                    
                                    // 
                                    label cellJ = dil[realCellI][interactingCells];
                                    
                                    if (0) {
                                    
                                    const scalar dh = pow(this->mesh().V()[cellJ], 1.0/3.0);
                                    const scalar thickness = dh * 0.5;
                                    const vector xc = this->mesh().C()[cellJ];
                                    
                                    scalar dist = mag(xc - xp) - rp;
                                    scalar alpha = 0;
                                    if (dist > thickness) {
                                        alpha = 0;
                                    } else if (dist < -thickness) {
                                        alpha = 1;
                                    } else {
                                        alpha = 0.5 - dist/thickness*0.5;
                                    }
                                    frac[cellJ] += alpha;
                                    } else {
                                        scalar alpha = calcPartInCellFrac(cellJ, xp, rp);
                                        frac[cellJ] += alpha;
                                    }
                                    //frac[cellJ] = 1.0;
                                }
                                
                                // Cell A itself
                                frac[realCellI] = 1.0;

                            }
                        }
                    }
                    
                    ibil_.receiveReferredData(pBufs, startOfRequests);
                    
                    // real-referred
                    {
                        // referred interaction list (RIL)
                        const labelListList &ril = ibil_.ril();
                        
                        List<IDLList<parcelType> > &referredParticles = 
                            ibil_.referredParticles();
                        
                        //List<DynamicList<parcelType*> > &cellOccupancy = this->cellOccupancy();
                        
                        // loop all referred cells
                        forAll(ril, refCellI)
                        {
                            IDLList<parcelType> &refCellRefParticles = 
                                referredParticles[refCellI];
                            
                            const labelList &realCells = ril[refCellI];
                            
                            // loop referred particles in referred cell
                            forAllIter
                            (
                                typename IDLList<parcelType>, 
                                refCellRefParticles,
                                referredParcel
                            )
                            {
                                const parcelType &pref = referredParcel();
                                const vector xp = pref.position();
                                const scalar dp = pref.d();
                                const scalar rp = dp * 0.5;
                                
                                // loop all real cells around referred cell
                                forAll(realCells, realCellI)
                                {
                                    label cellJ = realCells[realCellI];
                                    
                                    scalar alpha = calcPartInCellFrac(cellJ, xp, rp);
                                    frac[cellJ] += alpha;
                                    
                                    //frac[cellJ] += 1.0;
                                }
                            }
                        }
                    }
                    
                    
                    frac.correctBoundaryConditions();
                    
                    return tmp_ibfrac;
                }

        
        // Cloud evolution functions
            
            //- Store the current cloud state
            void storeState()
            {
                cloudCopyPtr_.reset
                (
                    static_cast<IBPartCloud<CloudType>*>
                    (
                        clone(this->name() + "Copy").ptr()
                    )
                );
            }
            
            //- Reset the current cloud to the previously stored state
            void restoreState()
            {
                cloudReset(cloudCopyPtr_());
                cloudCopyPtr_.clear();
            }
            
            //- Evolve the cloud
            void evolve()
            {
                if (this->solution().canEvolve())
                {
                    //
                    typename parcelType::template
                        TrackingData<IBPartCloud<CloudType> > td(*this);
                    
                    this->solve(td);
                }
            }
            
        // IO
            //- Print cloud information
            void info() 
            {
                CloudType::info();
                
                Info<< "   IB cloud" 
                    << nl;
            }


    // Member Operators

        //void operator=(const IBPartCloud&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        //friend Istream& operator>>(Istream&, IBPartCloud&);
        //friend Ostream& operator<<(Ostream&, const IBPartCloud&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "IBPartCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
