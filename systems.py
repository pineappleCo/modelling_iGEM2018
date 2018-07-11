import rates

class Systems():

  def __init__(self, imm_copy, col_copy, imm_str, col_str, imm_rbs, col_rbs):
    self.imm_copy_num = rates.copy_num[imm_copy]
    self.col_copy_num = rates.copy_num[col_copy]

    self.imm_promoter_str = rates.anderson_str[imm_str]
    self.col_promoter_str = rates.anderson_str[col_str]

    self.imm_rbs_affinity = rates.rbs_affinity[imm_rbs]
    self.col_rbs_affinity = rates.rbs_affinity[col_rbs]

    self.imm_transcription_rate = rates.imm_transcription_rate
    self.imm_deg_rate = rates.imm_deg_rate
    self.imm_translation_rate = rates.imm_translation_rate

    self.mRNA_deg_rate = rates.mRNA_deg_rate

    self.col_transcription_rate = rates.col_transcription_rate
    self.col_deg_rate = rates.col_deg_rate
    self.col_translation_rate = rates.col_translation_rate

    self.imm_transcription_rate_maxicell = rates.imm_transcription_rate_maxicell
    self.imm_deg_rate_maxicell = rates.imm_deg_rate_maxicell

    self.col_deg_rate_maxicell = rates.col_deg_rate_maxicell

  def imm_only(self, y, t):
     imm_mRNA = y[0]
     imm = y[1]
     #the model equations
     d_imm_mRNA_dt = (self.imm_copy_num * self.imm_promoter_str * self.imm_transcription_rate) - (self.mRNA_deg_rate * imm_mRNA)
     d_imm_dt = (self.imm_rbs_affinity * self.imm_translation_rate * imm_mRNA) - (self.imm_deg_rate * imm)
     return [d_imm_mRNA_dt, d_imm_dt]

  def imm_and_col(self, y, t):
     imm_mRNA = y[0]
     imm = y[1]
     col_mRNA = y[2]
     col = y[3]
     #the model equations
     d_imm_mRNA_dt = (self.imm_copy_num * self.imm_promoter_str * self.imm_transcription_rate) - (self.mRNA_deg_rate * imm_mRNA)
     d_imm_dt = (self.imm_rbs_affinity * self.imm_translation_rate * imm_mRNA) - (self.imm_deg_rate * imm)
     d_col_mRNA_dt = (self.col_copy_num * self.col_promoter_str * self.col_transcription_rate) - (self.mRNA_deg_rate * col_mRNA)
     d_col_dt = (self.col_rbs_affinity * self.col_translation_rate * col_mRNA) - (self.col_deg_rate * col)
     return [d_imm_mRNA_dt, d_imm_dt, d_col_mRNA_dt, d_col_dt]

  def imm_and_col_maxicell(self, y, t):
     imm_mRNA = y[0]
     imm = y[1]
     col_mRNA = y[2]
     col = y[3]
     # the model equations
     d_imm_mRNA_dt = (self.imm_copy_num * self.imm_promoter_str * self.imm_transcription_rate_maxicell) - (self.mRNA_deg_rate * imm_mRNA)
     d_imm_dt = (self.imm_rbs_affinity * self.imm_translation_rate * imm_mRNA) - (self.imm_deg_rate_maxicell * imm)
     d_col_mRNA_dt = (self.col_copy_num * self.col_promoter_str * self.col_transcription_rate) - (self.mRNA_deg_rate * col_mRNA)
     d_col_dt = (self.col_rbs_affinity * self.col_translation_rate * col_mRNA) - (self.col_deg_rate_maxicell * col)
     return [d_imm_mRNA_dt, d_imm_dt, d_col_mRNA_dt, d_col_dt]

  def imm_degrading_and_col(self, y, t):
     imm_mRNA = y[0]
     imm = y[1]
     col_mRNA = y[2]
     col = y[3]
     # the model equations
     d_imm_mRNA_dt = (self.imm_copy_num * self.imm_promoter_str * self.imm_transcription_rate_maxicell) - (self.mRNA_deg_rate * imm_mRNA)
     d_imm_dt = (self.imm_rbs_affinity * self.imm_translation_rate * imm_mRNA) - (self.imm_deg_rate * imm)
     d_col_mRNA_dt = (self.col_copy_num * self.col_promoter_str * self.col_transcription_rate) - (self.mRNA_deg_rate * col_mRNA)
     d_col_dt = (self.col_rbs_affinity * self.col_translation_rate * col_mRNA) - (self.col_deg_rate * col)
     return [d_imm_mRNA_dt, d_imm_dt, d_col_mRNA_dt, d_col_dt]
