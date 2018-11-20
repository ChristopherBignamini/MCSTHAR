#include "../Include/HadronDataSetLoader.h"
#include "../Include/HadronDataSet.h"
#include "../../../Utilities/Include/HadronizationException.h"
#include "../../../Utilities/Include/computeWaveFunctionStrangeComponent.h"

HadronDataSetLoader::HadronDataSetLoader(const string& i_hadronDataSetFileName,
                                         const double i_lightHadronMaxMass)
                                        :m_hadronDataSetAvailability(false)
{
    // Parse hadron data file
    auto_ptr< ::HadronDataSet> hadronDataSetFile;
    try
    {
        // Parse data set file
        hadronDataSetFile = HadronDataSet_(i_hadronDataSetFileName,
                                           xml_schema::flags::dont_validate);
    }
    catch(const xml_schema::exception& exception)
    {
        // Error during setup file parsing, throw exception
        string errorMessage("Error during hadron data set file loading, unable to parse the specified setup file. ");
        errorMessage += exception.what();
        throw HadronizationException(errorMessage,__FUNCTION__,421);
    }
    
    // Check number of available hadron data elements
    // TODO: add cross-check for sequence lenght (also in other loaders)
    const unsigned int numberOfHadronDataElements(hadronDataSetFile->hadronDataList().length());
    if(numberOfHadronDataElements==0)
    {
        // No hadron data element provided, throw exception
        throw HadronizationException("Error during hadron data set file loading, empty data set provided",
                                     __FUNCTION__,422);
    }
    m_hadronDataSet.resize(numberOfHadronDataElements);
    
    int idCode(0);
    unsigned int absIdCode(0);
    unsigned int flavour(0);
    unsigned int spinMultiplicity(0);
    int electricCharge(0);
    int strangeCharge(0);
    int charmCharge(0);
    int bottomCharge(0);
    int baryonicCharge(0);
    int topCharge(0);
    unsigned int heavyFlavoredConstituentNumber(0);
    double mass(0.);
    
    // Load single hadron data elements
    unsigned int numberOfLoadedHadronDataElements(0);
    for(unsigned int hadronDataIndex=0;hadronDataIndex<numberOfHadronDataElements;++hadronDataIndex)
    {
        
        // Retrieve hadron data
        const HadronDataElement& hadronData(hadronDataSetFile->hadronDataList().hadronData().at(hadronDataIndex));
        
        spinMultiplicity = hadronData.hadronSpinMultiplicity();
        if(spinMultiplicity==0)
        {
            // Wrong spin multiplicity value provided, throw exception
            throw HadronizationException("Error during hadron data set file loading, zero spin multiplicity provided",
                                         __FUNCTION__,423);
        }

        mass = hadronData.hadronMass();
        if(mass<=0.)
        {
            // Wrong mass value provided, throw exception
            throw HadronizationException("Error during hadron data set file loading, non positive mass provided",
                                         __FUNCTION__,424);
        }
        
        idCode = hadronData.hadronID();
        flavour = hadronData.hadronFlavorComposition();
        
        // Check if the current element is a hadron candidate
        absIdCode = std::abs(idCode);
        if(absIdCode>=111)
        {
            
            // Exclude K_long and K_short mesons
            if(absIdCode==130 || absIdCode==310)
            {
                continue;
            }
            
            // Anomalous hadrons (lighter than pi0) exclusion
            if(mass<0.1)
            {
                continue;
            }
            
            try
            {
                computeHadronCharges(flavour,
                                     idCode,
                                     electricCharge,
                                     strangeCharge,
                                     charmCharge,
                                     bottomCharge,
                                     topCharge,
                                     heavyFlavoredConstituentNumber);
                
            }
            catch(HadronizationException& ex)
            {
                // TODO: move to static string and error code for exception throw and catch
                if(ex.getReturnValue()==425)
                {
                    continue;
                }
                else
                {
                    throw ex;
                }
            }
            
            // Hadrons with more than one heavy flavored costituents are excluded
            // from the hadron set
            if(heavyFlavoredConstituentNumber>1)// TODO: this check should be implemented in HadronSet
            {
                continue;
            }
            
            // Hadrons with a top quark costituent are excluded from the hadron set
            if(topCharge>0)
            {
                continue;
            }
            
            // Light flavour hadron with mass larger than 1.8 GeV are excluded from
            // the hadron set
            if(heavyFlavoredConstituentNumber==0)
            {
                if(mass>i_lightHadronMaxMass)
                {
                    continue;
                }
            }

            const bool isMixingPresent(hadronData.mesonMixingData().present());
            double waveFunctionStrangeComponent(0.);
            unsigned int numberOfStrangeQuarks(std::abs(strangeCharge));
            if(flavour<100)
            {
                // Mesonic hadron type
                baryonicCharge = 0;
                if(isMixingPresent)
                {
                    if(strangeCharge!=0)
                    {
                        // Mixing provided for flavored meson, throw exception
                        throw HadronizationException("Error during hadron data set file loading, mixing data provided for flavored meson",
                                                     __FUNCTION__,427);
                    }
                    numberOfStrangeQuarks = 2;
                    
                    // TODO: avoid the following duplication
                    MesonSeries mesonFamily(Octect);
                    if(hadronData.mesonMixingData().get().mesonFamily() == "octect")
                    {
                    }
                    else if(hadronData.mesonMixingData().get().mesonFamily() == "singlet")
                    {
                        mesonFamily = Singlet;
                    }
                    else
                    {
                        // Not identified meson mixing family, throw exception// TODO: family or series?
                        throw HadronizationException("Error during hadron data set file loading, not identified family in mixing data",
                                                     __FUNCTION__,429);
                    }
                    waveFunctionStrangeComponent = computeWaveFunctionStrangeComponent(mesonFamily,
                                                                                       hadronData.mesonMixingData().get().mixingAngle());
                }
            }
            else
            {   
                // Baryonic hadron type
                if(waveFunctionStrangeComponent)
                {
                    // Mixing provided for baryon, throw exception
                    throw HadronizationException("Error during hadron data set file loading, mixing data provided for baryon",
                                                 __FUNCTION__,426);
                }
                
                if(idCode>0)
                {
                    baryonicCharge = 1;
                }
                else
                {
                    baryonicCharge = -1;
                }
            }
            
            m_hadronDataSet[numberOfLoadedHadronDataElements] = HadronData(idCode,
                                                                           strangeCharge,
                                                                           charmCharge,
                                                                           bottomCharge,
                                                                           static_cast<double>(electricCharge),
                                                                           static_cast<double>(baryonicCharge),
                                                                           mass,
                                                                           spinMultiplicity,
                                                                           numberOfStrangeQuarks,
                                                                           waveFunctionStrangeComponent);
            ++numberOfLoadedHadronDataElements;
        }
    }
    
    // Resize hadron data vector if required
    m_hadronDataSet.resize(numberOfLoadedHadronDataElements);
    
    m_hadronDataSetAvailability = true;
    return;
}

HadronDataSetLoader::~HadronDataSetLoader(void)
{
}

const vector<HadronData>& HadronDataSetLoader::getHadronDataSet(void) const
{
    if(m_hadronDataSetAvailability)
    {
        return m_hadronDataSet;
    }
    else
    {
        throw HadronizationException("Error during hadron data set retrieval, data not available",
                                     __FUNCTION__,428);
    }
}

void HadronDataSetLoader::computeHadronCharges(const unsigned int i_flavor,
                                               const int i_pid,
                                               int& io_electricCharge,
                                               int& io_strangeCharge,
                                               int& io_charmCharge,
                                               int& io_bottomCharge,
                                               int& io_topCharge,
                                               unsigned int& io_heavyFlavoredConstituentNumber) const
{
    io_electricCharge = 0;
    io_strangeCharge = 0;
    io_charmCharge = 0;
    io_bottomCharge = 0;
    io_topCharge = 0;
    io_heavyFlavoredConstituentNumber = 0;
    // Check new hadron composition
    if(i_flavor>10 || i_flavor<=555)
    {
        int quarkFlavour[3];
        unsigned int quarkNumber = 2;
        // TODO: change the if condition with "absFlavor<=55"
        //       in the final version. The current condition
        //       is useful only for validation vs old code.
        //       Moreover a better check could be implemented
        //       to avoid the inclusion of non standard states
        if(i_flavor<100)
        {
            // Retrieve meson quark flavors
            quarkFlavour[0] = i_flavor/10;
            quarkFlavour[1] = i_flavor - (i_flavor/10)*10;
        }
        else if(i_flavor>100)
        {
            // Retrieve baryon quark flavors
            quarkFlavour[0] = i_flavor/100;
            quarkFlavour[1] = ((i_flavor/10)*10 - (i_flavor/100)*100)/10;
            quarkFlavour[2] = i_flavor - (i_flavor/10)*10;
            quarkNumber = 3;
        }
        else
        {
            // Not identified hadron flavor
            stringstream hadronFlavorStr;
            hadronFlavorStr<<i_flavor;
            string excMessage = "Error during hadron set loading, not allowed hadron flavor composition (" + hadronFlavorStr.str() + ")";
            throw HadronizationException(excMessage,__FUNCTION__,425);
        }
        
        // Compute abelian charges
        int chargeSign = 1;
        for(unsigned int quarkIndex=0;quarkIndex<quarkNumber;++quarkIndex)
        {
            if(quarkIndex==1 && quarkNumber==2)
            {
                chargeSign = -1;
            }
            
            if(quarkFlavour[quarkIndex] == 1)
            {
                io_electricCharge -= chargeSign;
            }
            else if(quarkFlavour[quarkIndex] == 2)
            {
                io_electricCharge += 2*chargeSign;
            }
            else if(quarkFlavour[quarkIndex] == 3)
            {
                io_electricCharge -= chargeSign;
                io_strangeCharge -= chargeSign;
            }
            else if(quarkFlavour[quarkIndex] == 4)
            {
                io_electricCharge += 2*chargeSign;
                io_charmCharge += chargeSign;
                ++io_heavyFlavoredConstituentNumber;
            }
            else if(quarkFlavour[quarkIndex] == 5)
            {
                io_electricCharge -= chargeSign;
                io_bottomCharge -= chargeSign;
                ++io_heavyFlavoredConstituentNumber;
            }
            else if(quarkFlavour[quarkIndex] == 6)
            {
                io_electricCharge += 2*chargeSign;
                io_topCharge += chargeSign;
                ++io_heavyFlavoredConstituentNumber;
            }
            else if((quarkFlavour[quarkIndex] != 1) && (quarkFlavour[quarkIndex] != 2))
            {
                stringstream quarkFlavorStr;
                quarkFlavorStr<<quarkFlavour[quarkIndex];
                string excMessage = "Error during hadron set loading, not identified quark flavor (" + quarkFlavorStr.str() + ")";
                throw HadronizationException(excMessage,__FUNCTION__,425);
            }
            
        }
        
        // Change charge sign for anti-baryons
        if(quarkNumber == 3 && i_pid<0)
        {
            io_electricCharge *= -1;
            io_strangeCharge *= -1;
            io_charmCharge *= -1;
            io_bottomCharge *= -1;
            io_topCharge *= -1;
        }
        
    }
    else
    {
        // Not identified hadron flavor
        stringstream hadronFlavorStr;
        hadronFlavorStr<<i_flavor;
        string excMessage = "Error during hadron set loading, not allowed hadron flavor composition (" + hadronFlavorStr.str() + ")";
        throw HadronizationException(excMessage,__FUNCTION__,425);
    }
    
    // Normalize electric charge
    io_electricCharge /=3;
    
    return;
}
