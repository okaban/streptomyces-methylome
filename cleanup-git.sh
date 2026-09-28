#!/usr/bin/env bash
# epi-trans / streptomyces-methylome リポジトリ整理 (v3)
# 通常実行: cd /Users/okaban/bioinfo/rna-seq && bash cleanup-git.sh   (未コミットありはスキップ＋内訳表示)
# 強制実行: FORCE=1 bash cleanup-git.sh                                (未コミットも捨てて全削除)
set -u
cd /Users/okaban/bioinfo/rna-seq || exit 1
main="$(git rev-parse --show-toplevel)"
FORCE="${FORCE:-0}"

echo "==== BEFORE ===="
echo "worktrees: $(git worktree list | wc -l | tr -d ' ')  /  claude branches: $(git branch | grep -c 'claude/')  /  .git: $(du -sh .git | cut -f1)"
echo

removed=0; skipped=0
declare -A kinds
sample=0
while IFS= read -r w; do
  [ "$w" = "$main" ] && continue
  case "$w" in *"/.claude/worktrees/"*) ;; *) continue ;; esac
  dirty="$(git -C "$w" status --porcelain 2>/dev/null)"
  if [ -n "$dirty" ] && [ "$FORCE" != "1" ]; then
    skipped=$((skipped+1))
    while IFS= read -r line; do st="${line:0:2}"; kinds["$st"]=$(( ${kinds["$st"]:-0} + 1 )); done <<< "$dirty"
    if [ "$sample" -lt 3 ]; then echo "SKIP $w"; echo "$dirty" | head -4 | sed 's/^/    /'; sample=$((sample+1)); fi
  else
    git worktree remove "$w" 2>/dev/null || git worktree remove --force "$w" 2>/dev/null
    removed=$((removed+1))
  fi
done < <(git worktree list --porcelain | awk '/^worktree /{print $2}')

echo
echo "removed: $removed   skipped: $skipped"
if [ "$skipped" -gt 0 ]; then
  echo "-- スキップ分の未コミット内訳（ステータス別ファイル数。?? = 未追跡, M = 変更, D = 削除）--"
  for k in "${!kinds[@]}"; do echo "   [$k] ${kinds[$k]}"; done
  echo ">> 中身が不要なら強制削除:  FORCE=1 bash cleanup-git.sh"
fi

git worktree prune
git branch --merged main | grep 'claude/' | tr -d ' ' | xargs -r -n1 git branch -d 2>/dev/null
git gc --prune=now >/dev/null 2>&1

echo
echo "==== AFTER ===="
echo "worktrees: $(git worktree list | wc -l | tr -d ' ')  /  claude branches: $(git branch | grep -c 'claude/')  /  .git: $(du -sh .git | cut -f1)"
