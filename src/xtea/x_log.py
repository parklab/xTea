
class XLog():

    def open_file(self, sf_log):
        with open(sf_log, "w") as fout:#everytime start a new file
            fout.write("[Log starts:]\n")

        f_log=open(sf_log, "a+")
        return f_log

    def append_to_file(self, f_log, s_content):
        f_log.write(s_content+"\n")
        f_log.flush()

    def close_file(self, f_log):
        f_log.close()
