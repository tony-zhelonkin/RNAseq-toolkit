# Migration Guide: Transitioning Existing Projects to New Branching Strategy

This guide helps you migrate your existing projects that use the RNAseq-toolkit submodule from the old `dev` branch to the new project-specific branch workflow.

## Background

**Old workflow:** All projects directly used and modified the `dev` branch.

**New workflow:** Each project has its own `dev-{ProjectName}` branch, which can:
- Develop features independently
- Optionally contribute back to the shared `dev` branch
- Pull improvements from `dev` when ready

## Migration Steps for Each Project

### Step 1: Identify Your Projects

List all projects currently using RNAseq-toolkit:
- GVDRP1_prj (Already migrated to `dev-GVDRP1`)
- Project2 (To be migrated to `dev-Project2`)
- Project3 (To be migrated to `dev-Project3`)
- etc.

### Step 2: For Each Project, Create a Project-Specific Branch

Navigate to the RNAseq-toolkit repository (not the submodule in your project):

```bash
# Clone or navigate to standalone RNAseq-toolkit repo
cd /path/to/RNAseq-toolkit

# Ensure you have latest dev
git fetch origin
git checkout dev
git pull origin dev

# Create project-specific branch from dev
git checkout -b dev-Project2
git push -u origin dev-Project2

# Repeat for each project:
git checkout dev
git checkout -b dev-Project3
git push -u origin dev-Project3
```

### Step 3: Update Each Project's Submodule Configuration

For each parent project, update the submodule to track its specific branch:

```bash
# Navigate to parent project
cd /path/to/Project2

# Update .gitmodules file
# Add the line: branch = dev-Project2
# Under [submodule "path/to/RNAseq-toolkit"]

# Example .gitmodules after update:
# [submodule "01_Scripts/RNAseq-toolkit"]
#     path = 01_Scripts/RNAseq-toolkit
#     url = git@github.com:tony-zhelonkin/RNAseq-toolkit.git
#     branch = dev-Project2
```

Edit your `.gitmodules` file or use this command:

```bash
# Navigate to parent project root
cd /path/to/Project2

# Edit .gitmodules manually, or use:
git config -f .gitmodules submodule.01_Scripts/RNAseq-toolkit.branch dev-Project2

# Navigate to submodule and switch to project branch
cd 01_Scripts/RNAseq-toolkit
git fetch origin
git checkout dev-Project2

# Return to parent and commit the change
cd ../..
git add .gitmodules 01_Scripts/RNAseq-toolkit
git commit -m "Migrate RNAseq-toolkit submodule to dev-Project2 branch"
git push
```

### Step 4: Verify the Migration

Test that the submodule is properly configured:

```bash
# In parent project root
git submodule update --remote

# Navigate to submodule
cd 01_Scripts/RNAseq-toolkit

# Check current branch
git branch --show-current
# Should output: dev-Project2 (or your project's branch name)

# Verify tracking
git status
# Should show: "Your branch is up to date with 'origin/dev-Project2'"
```

## Migration Checklist

For each project:

- [ ] Create project-specific branch in RNAseq-toolkit repo (`dev-ProjectName`)
- [ ] Push branch to GitHub
- [ ] Update parent project's `.gitmodules` to specify branch
- [ ] Switch submodule to project-specific branch
- [ ] Commit and push parent project changes
- [ ] Verify submodule tracks correct branch
- [ ] Document which features (if any) should be contributed back to `dev`

## Post-Migration Workflow

### Regular Development

```bash
# In your project's submodule
cd path/to/RNAseq-toolkit
git checkout dev-ProjectName

# Make changes, commit, push
git add .
git commit -m "Add feature for ProjectName"
git push origin dev-ProjectName

# In parent project, commit submodule update
cd ../..
git add path/to/RNAseq-toolkit
git commit -m "Update RNAseq-toolkit: Add feature X"
git push
```

### Contributing Features to Shared `dev`

When you develop something useful for other projects:

```bash
# Push your changes to your project branch
cd path/to/RNAseq-toolkit
git push origin dev-ProjectName

# Create Pull Request on GitHub:
# Base: dev
# Compare: dev-ProjectName
#
# After review and merge, other projects can pull your improvements
```

### Getting Updates from Other Projects

```bash
# Periodically merge shared dev into your project branch
cd path/to/RNAseq-toolkit
git checkout dev-ProjectName
git fetch origin
git merge origin/dev
git push origin dev-ProjectName

# Update parent project
cd ../..
git add path/to/RNAseq-toolkit
git commit -m "Sync RNAseq-toolkit with latest dev improvements"
git push
```

## Troubleshooting

### Submodule shows detached HEAD

```bash
cd path/to/RNAseq-toolkit
git checkout dev-ProjectName
git pull origin dev-ProjectName
```

### Submodule not updating to correct branch

```bash
# Verify .gitmodules configuration
cat .gitmodules

# Force submodule to update
git submodule update --remote --force

# Or manually:
cd path/to/RNAseq-toolkit
git fetch origin
git checkout dev-ProjectName
```

### Conflicts when merging dev into project branch

```bash
cd path/to/RNAseq-toolkit
git checkout dev-ProjectName
git fetch origin
git merge origin/dev

# Resolve conflicts
# Edit conflicted files
git add .
git commit -m "Merge dev into dev-ProjectName: resolve conflicts"
git push origin dev-ProjectName
```

## FAQ

**Q: Do I have to contribute my changes back to `dev`?**

A: No! Only contribute features that would be useful for other projects. Project-specific code stays in your project branch.

**Q: Can I keep my project branch permanently separate?**

A: Yes. You can develop independently and only merge from `dev` when you want updates. You never have to merge back.

**Q: What if I want to switch my project back to the old `dev` branch?**

A: Just update `.gitmodules` to `branch = dev`, run `git submodule update --remote`, and commit. But this defeats the purpose of the new workflow.

**Q: How do I know when `dev` has new features I should pull?**

A: Check the RNAseq-toolkit repository's `dev` branch commit history or pull requests. Coordinate with other users if needed.

## Status Tracking

| Project | Branch Name | Migrated? | Migration Date | Notes |
|---------|-------------|-----------|----------------|-------|
| GVDRP1_prj | dev-GVDRP1 | ✅ Yes | 2025-11-19 | Initial migration |
| Project2 | dev-Project2 | ⬜ No | - | Pending |
| Project3 | dev-Project3 | ⬜ No | - | Pending |
| Project4 | dev-Project4 | ⬜ No | - | Pending |

Update this table as you migrate each project.
