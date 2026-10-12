#!/bin/bash
# Shared helpers for reading raw MinKNOW output (sourced by 02_filter_reads.sh and
# count_raw_reads.sh). Written for macOS bash 3.2 and BSD awk as well as Linux.
#
# A barcode folder may hold several copies of the same chunk, e.g.
#   ..._7.fastq.gz   ..._7.fastq   ..._7_sanitized.fastq   ..._7_log.txt
# Every read must be counted once, so each chunk is reduced to ONE file:
#   1. the .fastq.gz (MinKNOW original), if it is complete;
#   2. otherwise whichever copy holds the most complete FASTQ records
#      (ties go to .fastq.gz, then .fastq, then _sanitized.fastq).
# A damaged file (cut off when live basecalling stopped) is read up to the damage.

# reads <TAB> bases <TAB> damaged(0/1) <TAB> plain_text_named_gz(0/1)
fq_check () {
  local f="$1" gzbad=0 plain=0
  if [ "$(head -c 2 "$f" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then
    gzip -t "$f" 2>/dev/null || gzbad=1
  else
    case "$f" in *.gz) plain=1 ;; esac          # named .gz but not compressed
  fi
  { gzip -dcf "$f" 2>/dev/null || true; } | awk -v gzbad="$gzbad" -v plain="$plain" '
    NR % 4 == 1 { h = $0; next }
    NR % 4 == 2 { s = $0; next }
    NR % 4 == 3 { p = $0; next }
    { if (substr(h,1,1) == "@" && substr(p,1,1) == "+" && length(s) > 0 && length(s) == length($0)) { r++; b += length(s) }
      else { bad = 1; exit } }
    END { if (NR % 4 != 0) bad = 1
          printf "%d\t%d\t%d\t%d\n", r, b, (bad || gzbad) ? 1 : 0, plain }'
}

# Stream the complete records of a FASTQ (gzipped or not), stopping at the first damaged one.
fq_complete () {
  { gzip -dcf "$1" 2>/dev/null || true; } | awk '
    NR % 4 == 1 { h = $0; next }
    NR % 4 == 2 { s = $0; next }
    NR % 4 == 3 { p = $0; next }
    { if (substr(h,1,1) == "@" && substr(p,1,1) == "+" && length(s) > 0 && length(s) == length($0))
        print h "\n" s "\n" p "\n" $0
      else exit }'
}

# Chunk name shared by all copies of one chunk.
chunk_key () {
  local b; b=$(basename "$1")
  b=${b%.gz}; b=${b%.fastq}; b=${b%.fq}; b=${b%_sanitized}
  printf "%s" "$b"
}

# One line per chunk in a folder:
#   chunk <TAB> chosen_file <TAB> reads <TAB> bases <TAB> damaged <TAB> plain_gz <TAB> copies <TAB> note
fq_pick () {
  local d="$1" f key rank
  local cur="" best="" br=0 bb=0 bd=0 bp=0 ncopy=0 note="" done_chunk=0
  local r b dmg pl
  finish () {
    [ -n "$cur" ] || return 0
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$cur" "$best" "$br" "$bb" "$bd" "$bp" "$ncopy" "${note:-intact}"
  }
  for f in "$d"/*.fastq.gz "$d"/*.fq.gz "$d"/*.fastq "$d"/*.fq; do
    [ -e "$f" ] || continue
    case "$f" in *.gz) rank=1 ;; *_sanitized.*) rank=3 ;; *) rank=2 ;; esac
    printf "%s\t%s\t%s\n" "$(chunk_key "$f")" "$rank" "$f"
  done | sort -t "$(printf '\t')" -k1,1 -k2,2n | {
    while IFS="$(printf '\t')" read -r key rank f; do
      if [ "$key" != "$cur" ]; then
        finish
        cur="$key"; best=""; br=0; bb=0; bd=0; bp=0; ncopy=0; note=""; done_chunk=0
      fi
      ncopy=$((ncopy + 1))
      [ "$done_chunk" = 1 ] && continue           # an intact first copy was already taken
      IFS="$(printf '\t')" read -r r b dmg pl < <(fq_check "$f")
      if [ -z "$best" ] || [ "$r" -gt "$br" ]; then
        best="$f"; br="$r"; bb="$b"; bd="$dmg"; bp="$pl"
      fi
      if [ "$ncopy" = 1 ] && [ "$dmg" = 0 ]; then
        done_chunk=1                              # highest-priority copy is complete: use it
      elif [ "$dmg" = 1 ]; then
        note="${note:+$note; }$(basename "$f") damaged ($r complete reads)"
      fi
    done
    finish
  }
}
