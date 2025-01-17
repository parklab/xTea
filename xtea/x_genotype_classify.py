import os
import sys
import pandas as pd
from deepforest import CascadeForestClassifier
from scipy.io import arff

class GntpClassifier_DF21():
    def __init__(self):
        self.n_feature = 15
        return
####
    ####gnerate the arff file for training data (with label)
    def gnrt_training_arff_from_xTEA_output(self, sf_00_list, sf_01_list, sf_11_list, sf_arff, b_balance=False):
        with open(sf_arff, "w") as fout_arff:
            fout_arff.write("@RELATION\tinsgntp\n\n")
            fout_arff.write("@ATTRIBUTE\tlclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tpolyA\tNUMERIC\n")#total polyA

            # fout_arff.write("@ATTRIBUTE\trpolyA\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tlcov\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trcov\tNUMERIC\n")

            fout_arff.write("@ATTRIBUTE\tclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tfullmap\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tclipratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawlclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawrclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscordant\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tconcordant\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tindel\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tclass\t{0,1,2}\n\n")
            fout_arff.write("@ATTRIBUTE\tclass\t{1,2}\n\n") #two category
            fout_arff.write("@DATA\n")
####
            # with open(sf_00_list) as fin_00:
            #     for line in fin_00:
            #         sf_rslt = line.rstrip()
            #         l_features = self.load_in_feature_from_xTEA_output(sf_rslt)
            #         for rcd in l_features:
            #             fout_arff.write(",".join(rcd) + "\n") ##

            n_11=0
            with open(sf_11_list) as fin_11:
                for line in fin_11:
                    sf_rslt = line.rstrip()
                    l_features = self.load_in_feature_from_xTEA_output(sf_rslt)
                    for rcd in l_features:
                        fout_arff.write(",".join(rcd) + "\n")
                    n_11+=1

            n_01=0
            with open(sf_01_list) as fin_01:
                for line in fin_01:
                    if n_01>n_11 and b_balance==True:
                        break
                    sf_rslt = line.rstrip()
                    l_features = self.load_in_feature_from_xTEA_output(sf_rslt)
                    for rcd in l_features:
                        fout_arff.write(",".join(rcd) + "\n")
                    n_01+=1

####
    # ####train the model
    # def train_model(self, sf_arff, sf_model, f_ratio=0.3):
    #     data = arff.loadarff(sf_arff)
    #     df = pd.DataFrame(data[0])
    #     xVar = df.iloc[:, :self.n_feature]
    #     yVar = df.iloc[:, self.n_feature]
    #     yVar = yVar.astype('int')
    #
    #     X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=f_ratio)
    #     clf = RandomForestClassifier(n_jobs=-1, random_state=0, n_estimators=20)
    #     # clf = svm.SVC(kernel='linear')
    #     # clf=svm.SVC(kernel='rbf') #Gaussian Kernel
    #     clf.fit(X_train, y_train)
    #     with open(sf_model, 'wb') as file:
    #         pickle.dump(clf, file)
    #     preds = clf.predict(X_test)
    #     accuracy = accuracy_score(y_test, preds)
    #     print('Mean accuracy score: {0}'.format(round(accuracy)))
    #     tab = pd.crosstab(y_test, preds, rownames=['Actual Result'], colnames=['Predicted Result'])
    #     print(tab)

    ####
    ####
    ####Given xTEA output, generate the prepared arff file
    def prepare_arff_from_xTEA_output(self, sf_xTEA, sf_arff):
        with open(sf_arff, "w") as fout_arff:
            fout_arff.write("@RELATION\tinsgntp\n\n")
            fout_arff.write("@ATTRIBUTE\tlclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tpolyA\tNUMERIC\n")  # total polyA

            # fout_arff.write("@ATTRIBUTE\trpolyA\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tlcov\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trcov\tNUMERIC\n")

            fout_arff.write("@ATTRIBUTE\tclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tfullmap\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tclipratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawlclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawrclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscordant\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tconcordant\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tindel\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tclass\t{0,1,2}\n\n")
            fout_arff.write("@ATTRIBUTE\tclass\t{1,2}\n\n")
            fout_arff.write("@DATA\n")

            b_train = False
            l_features = self.load_in_feature_from_xTEA_output(sf_xTEA, b_train)
            for rcd in l_features:
                fout_arff.write(",".join(rcd) + "\n")

        data = arff.loadarff(sf_arff)
        df = pd.DataFrame(data[0])
        xVar = df.iloc[:, :self.n_feature]
        #yVar = df.iloc[:, self.n_feature]
        #yVar = yVar.astype('int')

        X=xVar.to_numpy()
        # X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=0.999999)
        # return X_test
        return X
####

    ####clf is the trained model
    def predict_for_site(self, sf_model, sf_xTEA, sf_new):
        rf_model_df21 = CascadeForestClassifier()
        rf_model_df21.load(sf_model)

        sf_arff = sf_xTEA + ".arff"
        # site_features=self.prepare_arff_from_xTEA_output_two_category(sf_xTEA, sf_arff)

        site_features = self.prepare_arff_from_xTEA_output(sf_xTEA, sf_arff)
        preds=None
        if len(site_features)>0:
            preds = rf_model_df21.predict(site_features)
        with open(sf_xTEA) as fin_xTEA, open(sf_new, "w") as fout_new:
            if None is preds:
                return
            i_idx = 0
            for line in fin_xTEA:
                sinfo = line.rstrip()
                s_gntp = "0/0"
                if preds[i_idx] == 1:
                    s_gntp = "0/1"
                elif preds[i_idx] == 2:
                    s_gntp = "1/1"
                sinfo += ("\t" + s_gntp + "\n")
                fout_new.write(sinfo)
                i_idx += 1

    ####
    def load_in_feature_from_xTEA_output(self, sf_xtea, b_train=True):
        l_all_features = []
        if os.path.isfile(sf_xtea)==False:
            return l_all_features
        with open(sf_xtea) as fin_xtea:
            for line in fin_xtea:
                fields = line.split()
                l_features = self._parser_features(fields)
                if b_train == True:
                    s_gntp = fields[-1]
                    if s_gntp == "FP" or s_gntp == "0/0":
                        s_gntp = "0"
                    elif s_gntp == "0/1" or s_gntp == "1/0":
                        s_gntp = "1"
                    elif s_gntp == "1/1":
                        s_gntp = "2"
                    l_features.append(s_gntp)
                else:
                    l_features.append("1")
                l_feature2 = []
                for tmp_rcd in l_features:
                    l_feature2.append(str(tmp_rcd))
                l_all_features.append(l_feature2)
        return l_all_features

####

####
    ####
    # parse out the features
    def _parser_features(self, l_fields):  #
        l_features = []
        f_lcov = float(l_fields[11])
        if f_lcov == 0:
            f_lcov = 0.0000000001
        f_rcov = float(l_fields[12])
        if f_rcov == 0:
            f_rcov = 0.0000000001
        l_features.append(float(l_fields[5])/f_lcov)#lclip-algn-on-consensus
        l_features.append(float(l_fields[6])/f_rcov)#rclip-algn-on-consensus
        l_features.append(float(l_fields[7])/f_lcov)#ldisc-algn-on-consensus
        l_features.append(float(l_fields[8])/f_rcov)#rdisc-algn-on-consensus
        l_features.append(float(l_fields[9])/f_lcov + float(l_fields[10])/f_rcov)#polyA

        # l_features.append()#right-polyA
        l_features.append(l_fields[11])#left-local-coverage
        l_features.append(l_fields[12])#right-local-coverage

        l_features.append(float(l_fields[35])/(f_lcov+f_rcov))#effective clip
        l_features.append(float(l_fields[36])/(f_lcov+f_rcov))#effective fully mapped

        f_clip_ratio=0
        if float(l_fields[36]) + float(l_fields[35]) > 0:
            f_clip_ratio = float(l_fields[35]) / (float(l_fields[36]) + float(l_fields[35]))  # clip ration
        l_features.append(f_clip_ratio)

        f_disc_ratio=0
        if (float(l_fields[39]) + float(l_fields[40]))>0:
            f_disc_ratio = float(l_fields[39]) / (float(l_fields[39]) + float(l_fields[40]))  # disc ratio
        l_features.append(f_disc_ratio)

        l_features.append(float(l_fields[37])/f_lcov)#raw left-clip
        l_features.append(float(l_fields[38])/f_rcov)# raw right-clip
        l_features.append(float(l_fields[39])/(f_lcov+f_rcov))#discordant pairs
        l_features.append(float(l_fields[40])/(f_lcov+f_rcov))# conordant pairs
        #l_features.append(float(l_fields[41]) / (f_lcov + f_rcov))  # indels reads

        self.n_feature = len(l_features)
        return l_features
#pkl_filename="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/genotyping/training_set_SSC/Genotyping/trained_model_ssc_py2_random_forest_two_category.pkl"
############################################################################################################
####
    def predict_for_site_vcf(self, sf_model, sf_xTEA_vcf, sf_new):
        rf_model_df21 = CascadeForestClassifier()
        rf_model_df21.load(sf_model)

        sf_arff = sf_xTEA_vcf + ".arff"
        site_features = self.prepare_arff_from_xTEA_vcf_output(sf_xTEA_vcf, sf_arff)
        preds=None
        if len(site_features)>0:
            preds = rf_model_df21.predict(site_features)
        with open(sf_xTEA_vcf) as fin_xTEA, open(sf_new, "w") as fout_new:
            if None is preds:
                return
            i_idx = 0
            for line in fin_xTEA:
                if len(line) > 0 and line[0] == "#":
                    fout_new.write(line.rstrip() + "\n")
                    continue

                l_fields = line.rstrip().split()
                s_info=""
                for s_term in l_fields[:-1]:  # here l_fields[-1] is the genotype field
                    s_info+=(s_term + "\t")
                s_gntp = "0/0"
                if preds[i_idx] == 1:
                    s_gntp = "0/1"
                elif preds[i_idx] == 2:
                    s_gntp = "1/1"
                fout_new.write(s_info+s_gntp+"\n")
                i_idx += 1

    def prepare_arff_from_xTEA_vcf_output(self, sf_xTEA_vcf, sf_arff):
        with open(sf_arff, "w") as fout_arff:
            fout_arff.write("@RELATION\tinsgntp\n\n")
            fout_arff.write("@ATTRIBUTE\tlclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trclipcns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tldisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trdisccns\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tpolyA\tNUMERIC\n")  # total polyA

            # fout_arff.write("@ATTRIBUTE\trpolyA\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tlcov\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trcov\tNUMERIC\n")

            fout_arff.write("@ATTRIBUTE\tclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tfullmap\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tclipratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscratio\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawlclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\trawrclip\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tdiscordant\tNUMERIC\n")
            fout_arff.write("@ATTRIBUTE\tconcordant\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tindel\tNUMERIC\n")
            #fout_arff.write("@ATTRIBUTE\tclass\t{0,1,2}\n\n")
            fout_arff.write("@ATTRIBUTE\tclass\t{1,2}\n\n")
            fout_arff.write("@DATA\n")

            b_train = False
            l_features = self.load_in_feature_from_xTEA_vcf_output(sf_xTEA_vcf, b_train)
            for rcd in l_features:
                fout_arff.write(",".join(rcd) + "\n")

        data = arff.loadarff(sf_arff)
        df = pd.DataFrame(data[0])
        xVar = df.iloc[:, :self.n_feature]
        #yVar = df.iloc[:, self.n_feature]
        #yVar = yVar.astype('int')

        X=xVar.to_numpy()
        # X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=0.999999)
        # return X_test
        return X

####
    # 1       101212112       .       G       <INS:ME:LINE1>  .       PASS    SVTYPE=INS:ME:LINE1;SVLEN=0;END=101212112;TSD=NULL;TSDLEN=-1;SUBTYPE=orphan_or_sibling_transduct
    # ion;TD_SRC=7:-1-78765277;STRAND=+;VAF=0.6382978723404256;LCLIP=14;RCLIP=0;LDISC=4;RDISC=3;LPOLYA=3;RPOLYA=0;LRAWCLIP=16;RRAWCLIP=7;VAF_CLIP=5;VAF_FMAP=0;VAF_DISC=13;VAF
    # _CONCORDNT=34;LDRC=0;LDNRC=0;RDRC=0;RDNRC=0;LCOV=40.77;RCOV=54.36;LD_AKR_RC=0;LD_AKR_NRC=0;RD_AKR_RC=0;RD_AKR_NRC=0;LC_CLUSTER=-1:-1;RC_CLUSTER=-1:-1;LD_CLUSTER=-1:-1;R
    # D_CLUSTER=78765277:78765277;NINDEL=4;CLIP_LEN=23:25:29:33:30;INS_INV=Not-5prime-inversion;REF_REP=not_in_LINE1_copy;GENE_INFO=not_gene_region   GT      1/1
    def load_in_feature_from_xTEA_vcf_output(self, sf_xtea_vcf, b_train=True):
        l_all_features = []
        if os.path.isfile(sf_xtea_vcf)==False:
            return l_all_features
        with open(sf_xtea_vcf) as fin_xtea:
            for line in fin_xtea:
                if len(line)>0 and line[0]=="#":
                    continue
                fields = line.split()
                l_features = self._parser_features_vcf(fields)
                if b_train == True:
                    s_gntp = fields[-1]
                    if s_gntp == "FP" or s_gntp == "0/0":
                        s_gntp = "0"
                    elif s_gntp == "0/1" or s_gntp == "1/0":
                        s_gntp = "1"
                    elif s_gntp == "1/1":
                        s_gntp = "2"
                    l_features.append(s_gntp)
                else:
                    l_features.append("1")
                l_feature2 = []
                for tmp_rcd in l_features:
                    l_feature2.append(str(tmp_rcd))
                l_all_features.append(l_feature2)
        return l_all_features

    def _parser_features_vcf(self, l_ori_fields):####
        l_features = []
        s_info=l_ori_fields[-3]
        l_fields=s_info.split(";")
        f_lclip, f_rclip, f_ldisc, f_rdisc, f_lpolyA, f_rpolyA, f_lrawclip, f_rrawclip =0,0,0,0,0,0,0,0
        f_vaf_clip, f_vaf_fmap, f_vaf_disc, f_vaf_concod, f_lcov, f_rcov=0,0,0,0,0,0

        for s_term in l_fields:
            l_term=s_term.split("=")
            if l_term[0]=="LCLIP":
                f_lclip=float(l_term[1])
            elif l_term[0]=="RCLIP":
                f_rclip = float(l_term[1])
            elif l_term[0]=="LDISC":
                f_ldisc = float(l_term[1])
            elif l_term[0]=="RDISC":
                f_rdisc = float(l_term[1])
            elif l_term[0]=="LPOLYA":
                f_lpolyA = float(l_term[1])
            elif l_term[0]=="RPOLYA":
                f_rpolyA = float(l_term[1])
            elif l_term[0]=="LRAWCLIP":
                f_lrawclip=float(l_term[1])
            elif l_term[0]=="RRAWCLIP":
                f_rrawclip=float(l_term[1])
            elif l_term[0]=="VAF_CLIP":
                f_vaf_clip=float(l_term[1])
            elif l_term[0]=="VAF_FMAP":
                f_vaf_fmap=float(l_term[1])
            elif l_term[0]=="VAF_DISC":
                f_vaf_disc=float(l_term[1])
            elif l_term[0]=="VAF_CONCORDNT":
                f_vaf_concod=float(l_term[1])
            elif l_term[0]=="LCOV":
                f_lcov = float(l_term[1])
            elif l_term[0]=="RCOV":
                f_rcov = float(l_term[1])


        if f_lcov == 0:
            f_lcov = 0.0000000001
        if f_rcov == 0:
            f_rcov = 0.0000000001

        l_features.append(f_lclip/f_lcov)
        l_features.append(f_rclip / f_rcov)
        l_features.append(f_ldisc / f_lcov)
        l_features.append(f_rdisc / f_rcov)
        l_features.append(f_lpolyA / f_lcov + f_rpolyA / f_rcov)
        l_features.append(f_lcov)
        l_features.append(f_rcov)
        l_features.append(f_vaf_clip/(f_lcov+f_rcov))
        l_features.append(f_vaf_fmap/(f_lcov+f_rcov))

        f_clip_ratio = 0
        if f_vaf_clip + f_vaf_fmap > 0:
            f_clip_ratio = f_vaf_clip / (f_vaf_clip + f_vaf_fmap)  # clip ration
        l_features.append(f_clip_ratio)

        f_disc_ratio = 0
        if f_vaf_disc + f_vaf_concod > 0:
            f_disc_ratio = f_vaf_disc / (f_vaf_disc + f_vaf_concod)  # disc ratio
        l_features.append(f_disc_ratio)

        l_features.append(f_lrawclip/f_lcov)
        l_features.append(f_rrawclip / f_rcov)
        l_features.append(f_vaf_disc / (f_lcov+f_rcov))
        l_features.append(f_vaf_concod / (f_lcov + f_rcov))

        self.n_feature = len(l_features)
        return l_features
