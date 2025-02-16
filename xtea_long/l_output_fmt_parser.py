import global_values

####this class saves some shared functions and names used in both 'LTransduction()' and 'LRepMasker()'
class LRegionNameBasicInfo():
    ####parse the region length
    def parse_region_fields(self, s_region):
        l_fields = s_region.split(":")
        i_start = int(l_fields[0])
        i_end = int(l_fields[1])
        i_len = i_end - i_start + 1
        return i_len

    def get_ins_info_from_qname(self, query_name):
        name_fields = query_name.split(global_values.SEPERATOR)
        ins_chrm = name_fields[0]
        ins_pos = name_fields[1]
        return ins_chrm, ins_pos

    def get_fields_from_str(self, s_region):
        l_fields = s_region.split(":")
        i_start = int(l_fields[0])
        i_end = int(l_fields[1])
        s_dir=l_fields[2]
        return i_start, i_end, s_dir

    def get_s_3mer(self):
        return "3prime"
    def get_s_5mer(self):
        return "5prime"
    def get_s_both_side(self):
        return "both"
    def get_s_NONE(self):
        return "None"
    def is_s_NONE(self, s_chk):
        return (s_chk == "None")
    #here need to add the set functions

class LInternalStructure():
    ####given two regions, return the internal structure
    def get_internal_structure(self, s_rg1, s_rg2, islack):
        b_inv = False
        b_del = False

        lbinfo = LRegionNameBasicInfo()
        if (lbinfo.is_s_NONE(s_rg1) == False) and (lbinfo.is_s_NONE(s_rg2) == False):
            (start1, end1, s_dir1) = lbinfo.get_fields_from_str(s_rg1)
            (start2, end2, s_dir2) = lbinfo.get_fields_from_str(s_rg2)

            if s_dir1 is not s_dir2:
                b_inv=True
            if (start2-end1) > islack:
                b_del=True
####
        # elif (s_rg1 is not None) and (s_rg2 is None):
        #     (start1, end1, s_dir1) = self.lbinfo.get_fields_from_str(s_rg1)
        # elif (s_rg1 is None) and (s_rg2 is not None):
        #     (start2, end2, s_dir2) = self.lbinfo.get_fields_from_str(s_rg2)
        s_structure = lbinfo.get_s_NONE()
        if b_inv==True and b_del==True:
            s_structure="internal_inversion_deletion"
        elif b_inv==True:
            s_structure = "internal_inversion"
        elif b_del==True:
            s_structure = "internal_deletion"
        return s_structure
####