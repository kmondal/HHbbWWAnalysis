import os
import sys
from copy import copy
from bamboo import treefunctions as op

#===============================================================================================#
#                                  SelectionObject class                                        #
#===============================================================================================#
class SelectionObject:
    """
        Helper class to hold a selection including :
        sel         : refine object
        selName     : name of the selection used in the refine
        yieldTitle  : Title for the yield table in LateX
        NB : copy() might be useful to be able to branch out the selection
    """
    def __init__(self,sel,selName,yieldTitle):
        self.sel        = sel
        self.selName    = selName
        self.yieldTitle = yieldTitle
    
    def refine(self,cut=None,weight=None,autoSyst=True):
        """
        Overload the refine function to avoid repeating the selName arg 
        """
        self.sel = self.sel.refine(name    = self.selName,
                                   cut     = cut,
                                   weight  = weight,
                                   autoSyst= autoSyst)
    def makeYield(self,yieldObject):
        """
        Record an entry in the yield table 
        """
        yieldObject.addYield(sel    = self.sel,
                             name   = self.selName,
                             title  = self.yieldTitle)

#===============================================================================================#
#                                        Selections                                             #
#===============================================================================================#

def makeSingleLeptonSelection(self,baseSel,plot_yield=False): 
    """
    Produces the requested lepton selection (encapsulated in SelectionObject class objects)
    Will produce a dict :
        - key = level required (Preselected, Fakeable, Tight and/or FakeExtrapolation)
        - value = list of SelectionObject class objects per channel [ElEl,MuMu,ElMu] 
    We start by leptons so no need to pass selObject
    Code organized such that selections are not repeated and hopefully optimzed the RooDataFrame
    """
    # Select level #
    sel_level = []   # What lepton type selection to perform (e.g. Tight needs Fakeable and Preselected)
    lepton_level = ["Preselected","Fakeable","Tight","FakeExtrapolation"]
    # What lepton type selection to save for plotting
    save_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in lepton_level and boolean]  
        # Make list of requested lepton selection in the arguments
    if len(save_level) == 0: save_level = lepton_level # If nothing asked, all will be done

    if "Tight" in save_level:
        sel_level.extend(["Preselected","Fakeable","Tight"]) # Tight needs Fakeable and Preselected
    if "FakeExtrapolation" in save_level:
        sel_level.extend(["Preselected","Fakeable","FakeExtrapolation"]) # FakeExtrapolation needs Fakeable and Preselected
    if "Fakeable" in save_level:
        sel_level.extend(["Preselected","Fakeable"]) # Fakeable needs Preselected
    if "Preselected" in save_level:
        sel_level.extend(["Preselected"]) 
            # TODO : find cleaner-> if args.Tight And args.FakeExtrapolation : will be redundancies (not a problem though)
    #--- Lambdas ---#
    # Fakeable lambdas #
    lambdaLowPtCut = lambda lepColl: lepColl[0].pt > 15 # leading above 15 GeV
    lambdaLeadingPtCut = lambda lepColl: lepColl[0].pt > 25 # leading above 25 GeV

    #Z-nominal mass
    Zmass = 91.1876
    # Mll cut lambdas #
    lambda_lowMllCut = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))
    # If any dilepton preselected below 12GeV, returns False
    lambda_outZ      = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(Zmass - 10., op.invariant_mass(dilep[0].p4, dilep[1].p4), Zmass + 10.)))
    # If any dilepton preselected within Z peak, returns False
    lambda_inZ       = lambda dilep : op.in_range(Zmass - 10., op.invariant_mass(dilep[0].p4, dilep[1].p4), Zmass + 10.)


    # Loose SF #
    ElLooseSF = lambda lepColl : [self.ttH_singleElectron_trigSF(lepColl[0])] + self.lambda_ElectronLooseSF(lepColl[0]) if self.is_MC else None
    MuLooseSF = lambda lepColl : [self.ttH_singleMuon_trigSF(lepColl[0])] + self.lambda_MuonLooseSF(lepColl[0]) if self.is_MC else None
#    ElLooseSF = lambda lepColl : self.lambda_ElectronLooseSF(lepColl[0]) if self.is_MC else None
#    MuLooseSF = lambda lepColl : self.lambda_MuonLooseSF(lepColl[0]) if self.is_MC else None    
    # Tight SF #
    ElTightSF = lambda lepColl : self.lambda_ElectronTightSF(lepColl[0]) if self.is_MC else None
    MuTightSF = lambda lepColl : self.lambda_MuonTightSF(lepColl[0]) if self.is_MC else None
    
    #--- Preselection ---#
    selectionDict = {}

    if "Preselected" in sel_level:
        ElPreSelObject = SelectionObject(sel          = baseSel,
                                         selName      = "HasElPreselected",
                                         yieldTitle   = "Preselected lepton (channel $e^+/e^-$)")
        MuPreSelObject = SelectionObject(sel          = baseSel,
                                         selName      = "HasMuPreselected",
                                         yieldTitle   = "Preselected lepton (channel $\mu^+/\mu^-$)")

        # PreSelectionObject # [PreSelIncl ????????]
        ElPreSelObject.refine(cut = [op.rng_len(self.electronsPreSel) >= 1])
        MuPreSelObject.refine(cut = [op.rng_len(self.muonsPreSel) >= 1])
        # Yield #
        if plot_yield:
            ElPreSelObject.makeYield(self.yieldPlots)
            MuPreSelObject.makeYield(self.yieldPlots)
        # Record if actual argument #
        if "Preselected" in save_level:
            selectionDict["Preselected"] = [ElPreSelObject,MuPreSelObject]

        #--- Fakeable ---#
        if "Fakeable" in sel_level:
            ElFakeSelObject = SelectionObject(sel         = ElPreSelObject.sel,
                                              selName     = "HasElFakeable",
                                              yieldTitle  = "Fakeable lepton (channel $e^+/e^-$)")
            MuFakeSelObject = SelectionObject(sel         = MuPreSelObject.sel,
                                              selName     = "HasMuFakeable",
                                              yieldTitle  = "Fakeable lepton (channel $\mu^+/\mu^-$)")

            # Selection : at least one fakeable dilepton #

####Commented out (Originally for DL)####
#            # rng_len needs to be done on unsliced container to make sure it is not empty for sliced container afterwards
#            ElElFakeSelObject.refine(cut = [op.rng_len(self.ElElDileptonFakeSel)>=1,
#                                            lambdaOSDilepton(self.ElElDileptonFakeSel[0]),
#                                            lambdaLowPtCutElEl(self.ElElDileptonFakeSel[0]),
#                                            lambdaLeadingPtCutElEl(self.ElElDileptonFakeSel[0])],
#                                     weight = ElElLooseSF(self.ElElDileptonFakeSel[0])) 
#            MuMuFakeSelObject.refine(cut = [op.rng_len(self.MuMuDileptonFakeSel)>=1,
#                                            lambdaOSDilepton(self.MuMuDileptonFakeSel[0]),
#                                            lambdaLowPtCutMuMu(self.MuMuDileptonFakeSel[0]),
#                                            lambdaLeadingPtCutMuMu(self.MuMuDileptonFakeSel[0])],
#                                     weight = MuMuLooseSF(self.MuMuDileptonFakeSel[0]))
#            ElMuFakeSelObject.refine(cut = [op.rng_len(self.ElMuDileptonFakeSel)>=1,
#                                            lambdaOSDilepton(self.ElMuDileptonFakeSel[0]),
#                                            lambdaLowPtCutElMu(self.ElMuDileptonFakeSel[0]),
#                                            lambdaLeadingPtCutElMu(self.ElMuDileptonFakeSel[0])],
#                                     weight = ElMuLooseSF(self.ElMuDileptonFakeSel[0]))
##            ElElFakeSelObject.refine(cut = [op.rng_len(self.ElElDileptonFakeSel)>=1,
##                                            lambdaOSDilepton(self.ElElDileptonFakeSel[0]),
##                                            lambdaLowPtCutElEl(self.ElElDileptonFakeSel[0]),
##                                            lambdaLeadingPtCutElEl(self.ElElDileptonFakeSel[0])])
##            MuMuFakeSelObject.refine(cut = [op.rng_len(self.MuMuDileptonFakeSel)>=1,
##                                            lambdaOSDilepton(self.MuMuDileptonFakeSel[0]),
##                                            lambdaLowPtCutMuMu(self.MuMuDileptonFakeSel[0]),
##                                            lambdaLeadingPtCutMuMu(self.MuMuDileptonFakeSel[0])])
##            ElMuFakeSelObject.refine(cut = [op.rng_len(self.ElMuDileptonFakeSel)>=1,
##                                            lambdaOSDilepton(self.ElMuDileptonFakeSel[0]),
##                                            lambdaLowPtCutElMu(self.ElMuDileptonFakeSel[0]),
##                                            lambdaLeadingPtCutElMu(self.ElMuDileptonFakeSel[0])])
##            ElEl_wSF = ElElLooseSF(self.ElElDileptonFakeSel[0])
##            MuMu_wSF = MuMuLooseSF(self.MuMuDileptonFakeSel[0])
##            ElMu_wSF = ElMuLooseSF(self.ElMuDileptonFakeSel[0]) 
##
##            from bamboo.analysisutils import forceDefine 
##            forceDefine(ElEl_wSF,ElElFakeSelObject.sel)
##            forceDefine(MuMu_wSF,MuMuFakeSelObject.sel)
##            forceDefine(ElMu_wSF,ElMuFakeSelObject.sel)
##
##            ElElFakeSelObject.selName += 'wgt'
##            MuMuFakeSelObject.selName += 'wgt'
##            ElMuFakeSelObject.selName += 'wgt'
##
##            ElElFakeSelObject.refine(weight = ElEl_wSF)
##            MuMuFakeSelObject.refine(weight = MuMu_wSF)
##            ElMuFakeSelObject.refine(weight = ElMu_wSF)
#############################################################
            ElFakeSelObject.refine(cut    = [op.rng_len(self.leadElectronFakeSel) >= 1, 
                                             lambdaLeadingPtCut(self.leadElectronFakeSel)],
                                   weight = ElLooseSF(self.leadElectronFakeSel))
            MuFakeSelObject.refine(cut    = [op.rng_len(self.leadMuonFakeSel) >= 1,
                                             lambdaLeadingPtCut(self.leadMuonFakeSel)],
                                   weight = MuLooseSF(self.leadMuonFakeSel))
            
            # Yield #
            if plot_yield:
                ElFakeSelObject.makeYield(self.yieldPlots)
                MuFakeSelObject.makeYield(self.yieldPlots)


            # Selection : Mll cut + Z peak exclusion in **preselected** leptons # 
            ElFakeSelObject.selName += "PreMllCut"
            MuFakeSelObject.selName += "PreMllCut"
           
            ElFakeSelObject.yieldTitle += " + $M_{ll}>12$"
            MuFakeSelObject.yieldTitle += " + $M_{ll}>12$"
            
            mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel),lambda_lowMllCut(self.ElMuDileptonPreSel)]

            ElEl_MllCut = copy(mllCut)
            MuMu_MllCut = copy(mllCut)
            ElMu_MllCut = copy(mllCut)

            if not self.args.NoZVeto: # No Mz cut in OS leptons
                # All the preselected leptons outside Z peak region within 10 GeV
                ElFakeSelObject.selName += "OutZ"
                MuFakeSelObject.selName += "OutZ"

                ElFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"
                MuFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"

                outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]

                ElEl_MllCut.extend(outZCut)
                MuMu_MllCut.extend(outZCut)
                ElMu_MllCut.extend(outZCut)

            ElFakeSelObject.refine(cut = ElEl_MllCut)
            MuFakeSelObject.refine(cut = MuMu_MllCut)

            # Yield #
            if plot_yield:
                ElFakeSelObject.makeYield(self.yieldPlots)
                MuFakeSelObject.makeYield(self.yieldPlots)


#            print("tauSel Size: %s   cleanedTauSel Size: %s"%(op.rng_len(self.tauSel),op.rng_len(self.tauCleanSel)))

                
            if not self.args.NoTauVeto: # nTau (isolated from fakeable leptons) = 0
                # All the preselected leptons outside Z peak region within 10 GeV
                ElFakeSelObject.selName += "noTau"
                MuFakeSelObject.selName += "noTau"

                ElFakeSelObject.yieldTitle += " + Tau Veto"
                MuFakeSelObject.yieldTitle += " + Tau Veto"

                ElFakeSelObject.refine(cut = [op.rng_len(self.tauCleanSel) == 0])
                MuFakeSelObject.refine(cut = [op.rng_len(self.tauCleanSel) == 0])
                
                # Yield #
                if plot_yield:
                    ElFakeSelObject.makeYield(self.yieldPlots)
                    MuFakeSelObject.makeYield(self.yieldPlots)
            
            # Record if actual argument #
            if "Fakeable" in save_level:
                selectionDict["Fakeable"] = [ElFakeSelObject,MuFakeSelObject]
                
            #--- Tight ---#
            if "Tight" in sel_level:
                ElTightSelObject = SelectionObject(sel          = ElFakeSelObject.sel,
                                                   selName      = "HasElTight",
                                                   yieldTitle   = "Tight lepton (channel $e^+/e^-$)")
                MuTightSelObject = SelectionObject(sel          = MuFakeSelObject.sel,
                                                   selName      = "HasMuTight",
                                                   yieldTitle   = "Tight lepton (channel $\mu^+/\mu^-$)")
                
                # Selection : at least one tight dilepton #
                ElTightSelObject.refine(cut    = [op.rng_len(self.leadElectronTightSel) == 1],
                                        weight = ElTightSF(self.leadElectronTightSel))
                MuTightSelObject.refine(cut    = [op.rng_len(self.leadMuonTightSel) == 1],
                                        weight = MuTightSF(self.leadMuonTightSel))

                # Yield #
                if plot_yield:
                    ElTightSelObject.makeYield(self.yieldPlots)
                    MuTightSelObject.makeYield(self.yieldPlots)
                
                # Record if actual argument #
                if "Tight" in save_level:
                    selectionDict["Tight"] = [ElTightSelObject,MuTightSelObject]


            #--- Fake extrapolated ---#
            if "FakeExtrapolation" in sel_level:
                ElFakeExtrapolationSelObject = SelectionObject(sel          = ElFakeSelObject.sel,
                                                               selName      = "HasElFakeExtrapolation",
                                                               yieldTitle   = "Fake extrapolated lepton (channel $e^+/e^-$)")
                MuFakeExtrapolationSelObject = SelectionObject(sel          = MuFakeSelObject.sel,
                                                               selName      = "HasMuFakeExtrapolation",
                                                               yieldTitle   = "Fake extrapolated lepton (channel $\mu^+/\mu^-$)")
                
                # Selection : at least one tight dilepton #
                ElFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.leadElectronFakeExtrapolationSel) >= 1])
                MuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.leadMuonFakeExrapolationSel) >= 1])
                
                # Yield #
                if plot_yield:
                    ElFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                    MuFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                    
                # Record if actual argument #
                if "FakeExtrapolation" in save_level:
                    selectionDict["FakeExtrapolation"] = [ElFakeExtrapolationSelObject,MuFakeExtrapolationSelObject]
    # Return # 
    return selectionDict

def makeAk4JetSelection(self,selObject,nJet,copy_sel=False,plot_yield=False):
    """
    Produces the Ak4 jet selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : at least two Ak4 jets
    """

    if copy_sel:
        selObject = copy(selObject)
        #    selObject.selName += "TwoAk4Jets"
        #    selObject.yieldTitle += " + Ak4 Jets $\geq 2$"
        #    selObject.refine(cut=[op.rng_len(self.ak4Jets)>=2])
    if nJet == 3: 
        selObject.selName += "Ak4JetsLoose"
        selObject.yieldTitle += " + Ak4 Jets Loose $\geq 2$"
    elif nJet == 4:
        selObject.selName += "Ak4JetsTight"
        selObject.yieldTitle += " + Ak4 Jets Tight $\geq 2$"
    selObject.refine(cut=[op.rng_len(self.ak4Jets) >= nJet])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel :
        return selObject 


def makeAk8JetSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the Ak8 jet selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : at least one Ak8 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "Ak8Jets"
    selObject.yieldTitle += " + Ak8 Jets $\geq 1$"
    selObject.refine(cut=[op.rng_len(self.ak8Jets) >= 1])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject


def makeExclusiveLooseResolvedJetComboSelection(self,selObject,nbJet,copy_sel=False,plot_yield=False):
    if copy_sel:
        selObject = copy(selObject)

    if nbJet == 0:
        AppliedSF = None
        '''
        #----- DY estimation from data -----#
        if self.args.DYDataEstimation1Btag or self.args.DYDataEstimation2Btag:
            if self.is_MC:
                print ("Warning : the DY reweighting is not applied on MC, will ignore that step")
            else:
                if "ElEl" in selObject.selName:
                    if self.DYReweightingElEl is None:
                        raise RuntimeError('DY reweighting for ElEl is not initialized')
                    AppliedDYReweighting = self.DYReweightingElEl(self.OSElElDileptonTightSel[0][0]) # Only on highest PT electron
                elif "MuMu" in selObject.selName:
                    if self.DYReweightingMuMu is None:
                        raise RuntimeError('DY reweighting for MuMu is not initialized')
                    AppliedDYReweighting = self.DYReweightingMuMu(self.OSMuMuDileptonTightSel[0][0]) # Only on highest PT muon
                else: # No DY reweighting in ElMu (marginal contribution)
                    AppliedDYReweighting = None
                AppliedSF = [AppliedDYReweighting] if AppliedDYReweighting is not None else None
                # weight = None works, but not weight = [None]
        '''
        selObject.selName += "ExclusiveResolved0b3j"
        selObject.yieldTitle += " + Exclusive Resolved (bjet = 0, ak4Jets >= 3)"
        selObject.refine(cut    = [op.rng_len(self.ak8BJets) == 0, 
                                   op.rng_len(self.ak4BJets) == 0,
                                   op.rng_len(self.ak4LightJets) >= 3],
                         weight = AppliedSF)

    elif nbJet == 1:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved1b2j"
        selObject.yieldTitle += " + Exclusive Resolved (bjet = 1, ak4Jets >= 2)"
        #exclusive 1b-Jet Selection (nAk4bJets=1 & nAk8bJets=0)
        #Applied the SFs as weight
        selObject.refine(cut    = [op.rng_len(self.ak8BJets) == 0, 
                                   op.rng_len(self.ak4BJets) == 1,
                                   op.rng_len(self.ak4LightJets) >= 2],
                         weight = AppliedSF)

    elif nbJet == 2:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved2b1j"
        selObject.yieldTitle += " + Exclusive Resolved (bjet = 2, ak4Jets >= 1)"
        '''        
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) >= 2,
                                   op.rng_len(self.ak8BJets) == 0],
                         weight = AppliedSF)
        
        if self.args.TTBarCR:
            selObject.selName += "MbbCut150"
            selObject.yieldTitle += " + $M_{bb}>150$"
            selObject.refine(cut = [op.invariant_mass(self.ak4BJets[0].p4,self.ak4BJets[1].p4)>150])
            if plot_yield:
                selObject.makeYield(self.yieldPlots)
        '''
        selObject.refine(cut    = [op.rng_len(self.ak8BJets) == 0,
                                   op.rng_len(self.ak4BJets) >= 2,
                                   op.rng_len(self.ak4LightJets) >= 1],
                         weight = None)
    else: raise RuntimeError ("Error in Jet Selection!!!")
    
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeExclusiveTightResolvedJetComboSelection(self,selObject,nbJet,copy_sel=False,plot_yield=False):
    if copy_sel:
        selObject = copy(selObject)

    if nbJet == 0:
        AppliedSF = None
        '''
        #----- DY estimation from data -----#
        if self.args.DYDataEstimation1Btag or self.args.DYDataEstimation2Btag:
            if self.is_MC:
                print ("Warning : the DY reweighting is not applied on MC, will ignore that step")
            else:
                if "ElEl" in selObject.selName:
                    if self.DYReweightingElEl is None:
                        raise RuntimeError('DY reweighting for ElEl is not initialized')
                    AppliedDYReweighting = self.DYReweightingElEl(self.OSElElDileptonTightSel[0][0]) # Only on highest PT electron
                elif "MuMu" in selObject.selName:
                    if self.DYReweightingMuMu is None:
                        raise RuntimeError('DY reweighting for MuMu is not initialized')
                    AppliedDYReweighting = self.DYReweightingMuMu(self.OSMuMuDileptonTightSel[0][0]) # Only on highest PT muon
                else: # No DY reweighting in ElMu (marginal contribution)
                    AppliedDYReweighting = None
                AppliedSF = [AppliedDYReweighting] if AppliedDYReweighting is not None else None
                # weight = None works, but not weight = [None]
        '''
        selObject.selName += "ExclusiveResolved0b4j"
        selObject.yieldTitle += " + Exclusive Resolved (0 bjet 4orMore ak4Jets)"
        selObject.refine(cut    = [op.rng_len(self.ak8BJets) == 0, 
                                   op.rng_len(self.ak4BJets) == 0,
                                   op.rng_len(self.ak4LightJets)>=4],
                         weight = AppliedSF)

    elif nbJet == 1:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved1b3j"
        selObject.yieldTitle += " + Exclusive Resolved (1 bjet 3orMore ak4Jets)"
        #exclusive 1b-Jet Selection (nAk4bJets=1 & nAk8bJets=0)
        #Applied the SFs as weight
        selObject.refine(cut    = [op.rng_len(self.ak8BJets) == 0,
                                   op.rng_len(self.ak4BJets) == 1,
                                   op.rng_len(self.ak4LightJets) >= 3],
                         weight = AppliedSF)

    elif nbJet == 2:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved2b2j"
        selObject.yieldTitle += " + Exclusive Resolved (2 bjets 2orMore ak4Jets)"
        '''
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) >= 2,
                                   op.rng_len(self.ak8BJets) == 0],
                         weight = AppliedSF)
        
        if self.args.TTBarCR:
            selObject.selName += "MbbCut150"
            selObject.yieldTitle += " + $M_{bb}>150$"
            selObject.refine(cut = [op.invariant_mass(self.ak4BJets[0].p4,self.ak4BJets[1].p4)>150])
            if plot_yield:
                selObject.makeYield(self.yieldPlots)
        '''
        selObject.refine(cut   = [op.rng_len(self.ak4BJets) >= 2,
                                  op.rng_len(self.ak8BJets) == 0,
                                  op.rng_len(self.ak4LightJets) >= 2],
                        weight = AppliedSF)

    else: raise RuntimeError ("Error in Jet Selection!!!")

    if plot_yield:
        selObject.makeYield(self.yieldPlots)

    if copy_sel:
        return selObject


def makeSemiBoostedHbbSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the inclusive boosted selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : At least one btagged Ak8 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    AppliedSF = None # TODO: correct at v7
    selObject.selName += "HbbBoosted"
    selObject.yieldTitle += " + Hbb Boosted"
    selObject.refine(cut    = [op.rng_len(self.ak8BJets) >=1,
                               op.rng_len(self.ak4JetsCleanedFromAk8b) >=1],
                     weight = AppliedSF)

    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject


def makeInclusiveBoostedSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the inclusive boosted selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : At least one btagged Ak8 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    AppliedSF = None # TODO: correct at v7
    selObject.selName += "InclusiveBoosted"
    selObject.yieldTitle += " + Inclusive Boosted"
    selObject.refine(cut    = [op.rng_len(self.ak8nonBJets) >= 1, 
                               op.rng_len(self.ak8BJets) >= 1],
                     weight = AppliedSF)

    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

