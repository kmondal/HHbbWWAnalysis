from bamboo.analysismodules import DataDrivenContribution

class DataDrivenFake(DataDrivenContribution):
    def usesSample(self, sampleName, sampleConfig):
        return not ('type' in sampleConfig.keys() and sampleConfig["type"]=="signal") 
        # Not to be done on signal samples
    def replacesSample(self, sampleName, sampleConfig):
        return False # No sample to be replaced : new contribution
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        # Fake contribution = reweighted data - reweighted MC (in Fake CR)
        modCfg = dict(sampleConfig)
        if 'type' in sampleConfig.keys() and sampleConfig["type"]=="signal":
            return {}
        if sampleConfig["group"]=="data":
            # Data : add as MC with correct normalization so that lumi*Xsec/generated-events = 1
            modCfg.update({"type": "mc",
                           "generated-events": lumi,
                           "cross-section": 1.,
                           "group": "Fake"})
        else:
            # MC : Change sign of contribution so that in the end Fake = Data-MC
            modCfg.update({"branching-ratio" : -1.,
                           "group": "Fake"})

        return modCfg

class DataDrivenDY(DataDrivenContribution):
    def usesSample(self, sampleName, sampleConfig):
        #return sampleConfig['group'] != 'DY' # Use data and non DY MC 
        return sampleConfig['group'] == 'data'
    def replacesSample(self, sampleName, sampleConfig):
        return sampleConfig['group'] == 'DY' # Replace DY by data evaluation
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        modCfg = dict(sampleConfig)
        if 'type' in sampleConfig.keys() and sampleConfig["type"]=="signal":
            return {}
        if sampleConfig["group"] == "data":
            # Data : add as MC with correct normalization so that lumi*Xsec/generated-events = 1
            modCfg.update({"type": "mc",
                           "generated-events": lumi,
                           "cross-section": 1.,
                           "group": "DYEstimation"})
#        else:
#            if sampleConfig["group"] == "DY":
#                return {}
#            else:
#                modCfg.update({"branching-ratio" : -1.,
#                               "group": "DYEstimation"})
        else:
            return {}

        return modCfg


