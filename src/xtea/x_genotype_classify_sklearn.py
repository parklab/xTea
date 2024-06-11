##06/02/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

#train the classification model
#predict based on the trained model

#this is a stand alone module for model training/prediction
import sys
import os
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.model_selection import train_test_split
from scipy.io import arff
import pandas as pd
#import pickle
import numpy
#from sklearn.metrics import accuracy_score

from skl2onnx import to_onnx
import onnxruntime as rt
from skl2onnx.common.data_types import FloatTensorType


#from sklearn import svm ####

class GntpClassifier_sklearn():
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
            fout_arff.write("@ATTRIBUTE\tpolyA\tNUMERIC\n")  #total polyA

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


    # ####train the model
    def train_model(self, sf_arff, sf_model, f_ratio=0.3):
        data = arff.loadarff(sf_arff)
        df = pd.DataFrame(data[0])
        xVar = df.iloc[:, :self.n_feature]
        yVar = df.iloc[:, self.n_feature]
        yVar = yVar.astype('int')

        X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=f_ratio)
        clf = RandomForestClassifier(n_jobs=-1, random_state=0, n_estimators=20)
        # clf = svm.SVC(kernel='linear')
        # clf=svm.SVC(kernel='rbf') #Gaussian Kernel
        clf.fit(X_train, y_train)
        onx = to_onnx(clf, X_train.to_numpy()[:1].astype(numpy.float32))
        with open("genotype_model_6_7_2024.onnx", "wb") as f:
            f.write(onx.SerializeToString())

        options = ort.SessionOptions()
        options.intra_op_num_threads = 1
        options.inter_op_num_threads = 1

        #with open(sf_model, 'wb') as file:
        #    pickle.dump(clf, file)
        preds = clf.predict(X_test)
        sess = rt.InferenceSession("genotype_model_6_7_2024.onnx", providers=["CPUExecutionProvider"])
        input_name = sess.get_inputs()[0].name
        label_name = sess.get_outputs()[0].name
        pred_onx = sess.run([label_name], {input_name: X_test.to_numpy().astype(numpy.float32)})[0]

        accuracy = accuracy_score(y_test, preds)
        accuracy_0 = accuracy_score(y_test, pred_onx)
        print('Mean accuracy score: {0}'.format(accuracy))
        print('Mean accuracy score ONX: {0}'.format(accuracy_0))
        tab = pd.crosstab(y_test, preds, rownames=['Actual Result'], colnames=['Predicted Result'])
        print(tab)

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
        yVar = df.iloc[:, self.n_feature]
        yVar = yVar.astype('int')
        # X_train, X_test, y_train, y_test = train_test_split(xVar, yVar, test_size=0.999999)
        # return X_test
        return xVar
        
    ####clf is the trained model
    def predict_for_site_onnx(self, file_model, sf_xTEA, sf_new):
        sf_arff = sf_xTEA + ".arff"

        options = rt.SessionOptions()
        options.intra_op_num_threads = 1
        options.inter_op_num_threads = 1

        sess = rt.InferenceSession(file_model, providers=["CPUExecutionProvider"])
        input_name = sess.get_inputs()[0].name
        label_name = sess.get_outputs()[0].name
            
        site_features = self.prepare_arff_from_xTEA_output(sf_xTEA, sf_arff)
        preds=None

        if len(site_features)>0:
            preds = sess.run([label_name], {input_name: site_features.to_numpy().astype(numpy.float32)})[0]

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


###

    ####clf is the trained model
    def predict_for_site(self, file_model, sf_xTEA, sf_new):
        sf_arff = sf_xTEA + ".arff"
        # site_features=self.prepare_arff_from_xTEA_output_two_category(sf_xTEA, sf_arff)

        with open(file_model, "rb") as f:
            rf_model = pickle.load(f)
            
        site_features = self.prepare_arff_from_xTEA_output(sf_xTEA, sf_arff)
        preds=None
        if len(site_features)>0:
            preds = rf_model.predict(site_features)
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
#

    def save_model(self, clf, pkl_filename):
        # Save to file in the current working directory
        # pkl_filename = "pickle_model.pkl"
        with open(pkl_filename, 'wb') as file:
            pickle.dump(clf, file)

    #
    def load_model_from_file(self, pkl_filename):
        # Load from file
        pickle_model = None
        with open(pkl_filename, 'rb') as file:
            if (sys.version_info > (3, 0)):
                #this is python3
                pickle_model = pickle.load(file, encoding='latin1')
                #pickle_model = pickle.load(file)
            else:
                #python2
                pickle_model = pickle.load(file)
        return pickle_model

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
####
