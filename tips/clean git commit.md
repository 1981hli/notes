# 清空 git commit

思想就是: 创建孤立分支; 删除 master 分支; 孤立分支改名为 master 分支; 推送

```bash
git checkout --orphan latest_branch
git add -A
git commit -am "First Commit"
git branch -D master
git branch -m master
git push -f origin master
```