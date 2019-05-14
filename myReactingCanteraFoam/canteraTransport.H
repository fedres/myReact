#include <vector>
#include <iostream>

#include "basicMultiComponentMixture.H"
//#include "cantera/transport.h"
//#include "cantera/thermo.h"
//#include "cantera/base/ctml.h"
#include "transport.h"
#include "MultiTransport.h"// to avoid name type error
#include "MixTransport.h"//same as line above
#include "thermo.h"
#include "base/ctml.h"

#pragma once

template<typename T>
T& force(T&& t){
   return t;
}


class CanteraTransport {
	private:
		const speciesTable& species_;
		const PtrList<volScalarField>& Y_;
//		const volScalarField& rho_;
		const volScalarField& p_;
		const volScalarField& T_;
        

		PtrList<volScalarField> diCoeffs_;
		PtrList<volScalarField> dtiCoeffs_;
	    PtrList<volScalarField> partialEnthalpies_;
		std::unique_ptr<volScalarField> kappa_;
		std::unique_ptr<volScalarField> mu_;

		std::unique_ptr<Cantera::ThermoPhase> canteraGas_;//check deprecated error
		Cantera::MixTransport canteraTran_;             //check does not name a type
		Cantera::MultiTransport canteraMultiTran_;
		std::vector<int> canteraInds_;
		std::vector<doublereal> molarmasses_;
		std::vector<doublereal> temp_;
        
        template<class FT, class YT, class ListT>
//        void updateStateVariables(const FT& T, const FT& rho, 
//             YT& Y, ListT& diCoeffs, ListT& dtiCoeffs, ListT& partialEnthalpies, FT& lambda)
        void updateStateVariables(const FT& T, const FT& p, 
             YT& Y, ListT& diCoeffs, ListT& dtiCoeffs, ListT& partialEnthalpies, FT& kappa, FT& mu)

        {
            forAll(T, cellI) {
                //Update the cantera state
                canteraGas_->setState_TP(T[cellI], p[cellI]);

                forAll(Y_, i) {
                    temp_[canteraInds_[i]] = Y[i][cellI];
                }
                canteraGas_->setMassFractions(temp_.data());
                
                //Update the diffusion coefficients
                canteraTran_.getMixDiffCoeffsMass(temp_.data());

                forAll(diCoeffs, i) {
                    diCoeffs[i][cellI] = temp_[canteraInds_[i]];
                }

                //Update the partial enthalpies
                canteraGas_->getPartialMolarEnthalpies(temp_.data());
                forAll(partialEnthalpies, i) {
                    partialEnthalpies[i][cellI] = temp_[canteraInds_[i]]/molarmasses_[canteraInds_[i]];
                }

                //Update the thermal diffusion coefficients
                canteraMultiTran_.getThermalDiffCoeffs(temp_.data());
                forAll(dtiCoeffs, i) {
                    dtiCoeffs[i][cellI] = temp_[canteraInds_[i]];
                }

                //Update the mixture thermal conductivity
                kappa[cellI] = canteraTran_.thermalConductivity();

                //Update the mixture viscosity
                mu[cellI] = canteraTran_.viscosity();
            }
        }


	public:
		CanteraTransport( 
			basicMultiComponentMixture& mix,
			const fvMesh& mesh
			) : 
				species_(mix.species()),
				Y_(mix.Y()),
//				rho_(mesh.lookupObject<volScalarField>("rho")), 
				p_(mesh.lookupObject<volScalarField>("p")), 
				T_(mesh.lookupObject<volScalarField>("T")),
				canteraInds_(species_.size()),
				molarmasses_(species_.size()),
				temp_(species_.size())
		{
			//===Initialise the cantera classes===
			cout << "Creating CanteraTransport model" << std::endl;
			const IOdictionary& thermoProperties = 
				mesh.lookupObject<IOdictionary>("thermophysicalProperties");
			
			//Generate the .cti file from the chemkin files	(chemistry is included but
			//presently not used
			if (!thermoProperties.found("CTMLFile")) { 
				fileName ckFile = thermoProperties.lookup("CHEMKINFile");
				fileName ckThermoFile = thermoProperties.lookup("CHEMKINThermoFile");
				fileName ckTransportFile = thermoProperties.lookup("CHEMKINTransportFile");
				fileName ctiFile = fileName(ckFile).lessExt() + ".cti";
				
				Cantera::ck2cti(ckFile.expand(), ckThermoFile.expand(), ckTransportFile.expand());
				
				canteraGas_.reset(Cantera::newPhase(ctiFile.expand(), "gas"));
			} else {
				fileName ctmlFile = thermoProperties.lookup("CTMLFile");
				canteraGas_.reset(Cantera::newPhase(ctmlFile.expand(), "gas"));
			}
			
			canteraGas_->getMolecularWeights(molarmasses_.data());

			//===Mapping between species indices in OpenFOAM and Cantera
			for(int i=0; i<species_.size(); i++)
				canteraInds_[i] = canteraGas_->speciesIndex(species_[i]);

			canteraTran_.init(canteraGas_.get());
			canteraMultiTran_.init(canteraGas_.get());

			
			Switch writeCoeffs = thermoProperties.lookupOrDefault("writeCoeffs", false);

			//===Initialise the OpenFOAM fields

			forAll(Y_, i) {
				diCoeffs_.append(new volScalarField
				(
					IOobject
					(
						"D_"+species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						writeCoeffs ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("D", dimLength*dimLength/dimTime, 0.0)
				));

				dtiCoeffs_.append(new volScalarField
				(
					IOobject
					(
						"DT_"+species_[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						writeCoeffs ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("DT", dimMass/dimLength/dimTime, 0.0)
				));
				
				partialEnthalpies_.append(new volScalarField
				(
					IOobject
					(
						"H_"+species_[i],
						mesh.time().timeName(),
						mesh
					),
					mesh,
					dimensionedScalar("H", dimEnergy/dimMass, 0.0)
				));
			}
			
			kappa_.reset(new volScalarField
			(
				IOobject
				(
					"kappa",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					writeCoeffs ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("kappa", dimPower/dimLength/dimTemperature,0.0)
			));

			mu_.reset(new volScalarField
			(
				IOobject
				(
					"mu",
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					writeCoeffs ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("mu", dimPressure*dimTime,0.0)
			));

			updateState();
			
			if (writeCoeffs) 
                        {
				kappa_->write();
                                mu_->write();

				forAll(Y_, i) 
				{
				   diCoeffs_[i].write();
				   dtiCoeffs_[i].write();
				}
			}
		}

				
		void updateState()
		{   //Update the internal volume field
//            updateStateVariables(T_, rho_, Y_, diCoeffs_, dtiCoeffs_, partialEnthalpies_, *lambda_);
            updateStateVariables(T_, p_, Y_, diCoeffs_, dtiCoeffs_, partialEnthalpies_, *kappa_, *mu_);
            //Update the boundary fields
            const volScalarField::Boundary& Tbf = T_.boundaryField();

            forAll(Tbf, patchi) {
                label n = species_.size();
                //now get the faces for each boundary patch of each field
                const fvPatchScalarField& pTbf = Tbf[patchi];
//                const fvPatchScalarField& prhobf = rho_.boundaryField()[patchi];
                 const fvPatchScalarField& ppbf = p_.boundaryField()[patchi];
                       fvPatchScalarField& pkappabf = kappa_->boundaryFieldRef()[patchi];// fixed using non constant acess
                       fvPatchScalarField& pmubf =   mu_->boundaryFieldRef()[patchi];  //same as above


                UPtrList<const fvPatchScalarField> pYbf(n);
                UPtrList< fvPatchScalarField> pdiCoeffsbf(n);
                UPtrList< fvPatchScalarField> pdtiCoeffsbf(n);
                UPtrList< fvPatchScalarField> ppartialEnthalpiesbf(n);

                forAll(Y_, i) {
                    pYbf.set        (i, &(Y_[i].boundaryField()[patchi]));
                    pdiCoeffsbf.set (i, &(diCoeffs_[i].boundaryFieldRef()[patchi]));//error
                    pdtiCoeffsbf.set(i, &(dtiCoeffs_[i].boundaryFieldRef()[patchi]));//error
                    ppartialEnthalpiesbf.set(i, &(partialEnthalpies_[i].boundaryFieldRef()[patchi]));//error
                }
                
//                updateStateVariables(pTbf, prhobf, pYbf, pdiCoeffsbf, pdtiCoeffsbf, ppartialEnthalpiesbf, plambdabf);
                updateStateVariables(pTbf, ppbf, pYbf, pdiCoeffsbf, pdtiCoeffsbf, ppartialEnthalpiesbf, pkappabf, pmubf);
                
            }
            
        }

		const PtrList<volScalarField>& partialEnthalpies()
		{
			return partialEnthalpies_;
		}

		const PtrList<volScalarField>& diffCoeffs()
		{
			return diCoeffs_;
		}	

		const volScalarField& kappa()
		{
			return *kappa_;
		}

		const volScalarField& mu()
		{
			return *mu_;
		}
		
		const PtrList<volScalarField>& diffTCoeffs()
		{
			return dtiCoeffs_;
		}
        
};

