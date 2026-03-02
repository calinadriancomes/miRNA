# ===========================
# R/02_io.R
# Excel I/O
# ===========================

.safe_read_sheet <- function(path, sheet = NULL) {
  if (is.null(sheet)) {
    # default: first sheet
    readxl::read_excel(path, sheet = 1)
  } else {
    readxl::read_excel(path, sheet = sheet)
  }
}

load_all_data <- function(smcp_path, mirna_path) {
  # SMPC.xlsx: user said first sheet is SMPC, second is miRNA
  smcp_sheetnames <- readxl::excel_sheets(smcp_path)

  # heuristic: pick "SMPC" and "miRNA" if exist
  smcp_sheet <- if ("SMPC" %in% smcp_sheetnames) "SMPC" else smcp_sheetnames[1]
  miRNA_sheet_in_smcp <- if ("miRNA" %in% smcp_sheetnames) "miRNA" else {
    if (length(smcp_sheetnames) >= 2) smcp_sheetnames[2] else smcp_sheetnames[1]
  }

  smcp <- .safe_read_sheet(smcp_path, sheet = smcp_sheet)
  smcp_mirna <- .safe_read_sheet(smcp_path, sheet = miRNA_sheet_in_smcp)

  # miRNA_0.xlsx: typically one sheet
  mirna_sheetnames <- readxl::excel_sheets(mirna_path)
  mirna <- .safe_read_sheet(mirna_path, sheet = mirna_sheetnames[1])

  list(
    smcp = smcp,
    smcp_mirna = smcp_mirna,
    mirna = mirna
  )
}