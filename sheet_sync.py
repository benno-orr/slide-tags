"""
sheet_sync.py -- pull/push slide-tags.csv <-> Google Sheet

Usage:
    python3 sheet_sync.py pull   # sheet is source of truth; overwrite local CSV
    python3 sheet_sync.py push   # local CSV is source of truth; update sheet

Registry entry (for Claude sheet-sync agent):
    table:     slide-tags
    sheet_id:  1ctXseMOjzDodmE581xLaFzc_dJ5I4svlMGTofDjNRUU
    worksheet: Sheet1
    local:     /n/data1/hms/scrb/chen/lab/bco/pipelines/slide-tags/slide-tags.csv
    key:       composite (xid, x_spl)   -- xBO368 appears twice with x_spl=1 and x_spl=2
"""

import gspread, csv, shutil, sys
from datetime import datetime
from google.oauth2.service_account import Credentials

KEY   = '/home/beo703/.config/gspread/bco-gs-scripts-3146376045cb.json'
SHEET_ID  = '1ctXseMOjzDodmE581xLaFzc_dJ5I4svlMGTofDjNRUU'
WORKSHEET = 'Sheet1'
LOCAL     = '/n/data1/hms/scrb/chen/lab/bco/pipelines/slide-tags/slide-tags.csv'
KEY_COLS  = ('xid', 'x_spl')   # composite key


def backup(path):
    dst = f"{path}.bak-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    shutil.copy(path, dst)
    return dst


def connect():
    creds = Credentials.from_service_account_file(
        KEY, scopes=['https://www.googleapis.com/auth/spreadsheets'])
    gc = gspread.authorize(creds)
    ws = gc.open_by_key(SHEET_ID).worksheet(WORKSHEET)
    return ws


def trim_blank_rows(rows):
    while rows and all(c.strip() == '' for c in rows[-1]):
        rows.pop()
    return rows


def pull():
    ws = connect()
    rows = trim_blank_rows(ws.get_all_values())
    bkp = backup(LOCAL)
    with open(LOCAL, 'w', newline='') as f:
        csv.writer(f).writerows(rows)
    print(f"pull: wrote {len(rows)} rows ({len(rows[0])} cols)")
    print(f"backup: {bkp}")


def push():
    ws = connect()
    sheet_rows = trim_blank_rows(ws.get_all_values())
    with open(LOCAL) as f:
        csv_rows = list(csv.reader(f))

    sheet_header = sheet_rows[0]
    csv_header   = csv_rows[0]
    sheet_col = {h: i for i, h in enumerate(sheet_header) if h.strip()}
    csv_col   = {h: i for i, h in enumerate(csv_header)   if h.strip()}
    shared    = set(sheet_col) & set(csv_col)

    def make_key(row, col_map):
        return tuple((row[col_map[k]].strip() if k in col_map and col_map[k] < len(row) else '')
                     for k in KEY_COLS)

    sheet_by_key = {}
    for i, r in enumerate(sheet_rows[1:], start=2):
        k = make_key(r, sheet_col)
        if k[0]:   # xid must be non-empty
            sheet_by_key[k] = (i, r)

    updates, updated, appended = [], [], []
    for r in csv_rows[1:]:
        k = make_key(r, csv_col)
        if not k[0]:
            continue
        if k in sheet_by_key:
            row_num, srow = sheet_by_key[k]
            for col in shared:
                si, ci = sheet_col[col], csv_col[col]
                lv = r[ci].strip()    if ci < len(r)    else ''
                sv = srow[si].strip() if si < len(srow) else ''
                if lv != sv:
                    updates.append(gspread.Cell(row=row_num, col=si+1, value=lv))
            if any(u.row == row_num for u in updates):
                updated.append(k)
        else:
            new_row = [''] * len(sheet_header)
            for col in shared:
                si, ci = sheet_col[col], csv_col[col]
                new_row[si] = r[ci] if ci < len(r) else ''
            n = len(sheet_rows) + len(appended) + 1
            end = gspread.utils.rowcol_to_a1(1, len(sheet_header))[:-1]
            ws.update(values=[new_row], range_name=f'A{n}:{end}{n}',
                      value_input_option='USER_ENTERED')
            appended.append(k)

    if updates:
        ws.update_cells(updates, value_input_option='USER_ENTERED')
    print(f"push: updated={len(set(updated))} appended={len(appended)} unchanged={len(csv_rows)-1-len(set(updated))-len(appended)}")


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'pull'
    if cmd == 'pull':
        pull()
    elif cmd == 'push':
        push()
    else:
        print(f"Unknown command: {cmd}. Use 'pull' or 'push'.")
        sys.exit(1)
