# Bioconductor Submission Guide for DESandbox

## Pre-Submission Checklist

### ✅ Completed
- [x] Package passes R CMD check with 0 errors, 0 warnings
- [x] All functions documented with roxygen2
- [x] Unit tests included (25 tests passing)
- [x] Vignette included and builds successfully
- [x] Valid DESCRIPTION file with biocViews
- [x] MIT license with LICENSE file
- [x] Proper imports and dependencies declared
- [x] Examples provided for main functions

### 📋 Before Submission

1. **Create GitHub Repository**
   ```bash
   # Initialize git if not already done
   git init
   git add .
   git commit -m "Initial commit - DESandbox v0.1.0"
   
   # Create repo on GitHub and push
   git remote add origin https://github.com/YOUR_USERNAME/DESandbox.git
   git branch -M main
   git push -u origin main
   ```

2. **Update DESCRIPTION file**
   - Replace placeholder author information:
     ```r
     Authors@R: c(
         person("Your", "Name", 
                email = "your.email@example.com", 
                role = c("aut", "cre"),
                comment = c(ORCID = "0000-0000-0000-0000"))
     )
     ```
   - Add your actual name, email, and ORCID ID

3. **Review README.md**
   - Ensure installation instructions are clear
   - Add badges (optional but recommended)
   - Include quick start example

4. **Final Check**
   ```r
   # In R console
   devtools::check()
   BiocCheck::BiocCheck(".")  # Install BiocCheck if needed
   ```

## Submission Process

### Step 1: Register on Bioconductor

1. Go to https://git.bioconductor.org/
2. Create an account if you don't have one
3. Set up SSH keys for authentication

### Step 2: Submit via GitHub

Bioconductor uses a GitHub-based submission process:

1. **Open a new issue** at: https://github.com/Bioconductor/Contributions/issues/new

2. **Use this template**:
   ```
   Title: DESandbox: Unified Interface for Differential Expression Analysis
   
   Package Name: DESandbox
   Version: 0.1.0
   
   GitHub Repository: https://github.com/YOUR_USERNAME/DESandbox
   
   Maintainer: Your Name <your.email@example.com>
   
   Brief Description:
   DESandbox provides a unified interface for running and comparing differential 
   expression analyses using DESeq2, edgeR, and limma-voom. It standardizes outputs,
   enables cross-method comparison, and integrates gene set enrichment analysis.
   Includes publication-quality visualization tools with Nature Genetics styling.
   
   Additional Notes:
   - Package passes R CMD check with 0 errors, 0 warnings
   - Includes comprehensive vignette and 25 unit tests
   - All functions fully documented
   ```

### Step 3: Review Process

The Bioconductor team will:

1. **Automated checks**: BiocCheck will run automatically
2. **Manual review**: Reviewers will examine:
   - Code quality and style
   - Documentation completeness
   - Vignette quality
   - Test coverage
   - Dependency appropriateness

3. **Feedback**: Address reviewer comments by:
   - Pushing fixes to your GitHub repo
   - Commenting on the issue thread
   - The package bot will auto-update

### Step 4: Common Issues to Address

Based on your current package:

1. **biocViews**: Your current biocViews are good, but consider adding more specific ones:
   ```
   biocViews: GeneExpression, DifferentialExpression, RNASeq, Sequencing,
       Transcriptomics, Normalization, Visualization, MultipleComparison,
       QualityControl, ImmunoOncology
   ```

2. **Version Numbering**: Bioconductor uses odd minor versions for development:
   - Change from `0.1.0` to `0.99.0` (pre-release version)
   - After acceptance, it becomes `1.0.0` in the next release

3. **NEWS File**: Create a NEWS.md file:
   ```markdown
   # DESandbox 0.99.0
   
   ## New Features
   * Initial Bioconductor submission
   * Unified interface for DESeq2, edgeR, and limma-voom
   * Publication-quality visualizations
   * Cross-method comparison tools
   * Gene set enrichment integration
   ```

4. **man/DESandbox-package.Rd**: Create package-level documentation

## Quick Commands

```bash
# Update version to pre-release
# Edit DESCRIPTION: Version: 0.99.0

# Rebuild documentation
Rscript -e "roxygen2::roxygenise()"

# Final check
Rscript -e "devtools::check()"

# Install BiocCheck and run
Rscript -e "BiocManager::install('BiocCheck')"
Rscript -e "BiocCheck::BiocCheck('.')"

# Build tarball for local testing
R CMD build .
R CMD check DESandbox_0.99.0.tar.gz --as-cran
```

## Timeline

- **Submission**: Immediate (after completing checklist)
- **Initial bot checks**: Within hours
- **Reviewer assignment**: 1-2 weeks
- **Review process**: 2-6 weeks (depends on feedback cycles)
- **Acceptance**: Added to next Bioconductor release (April or October)

## Resources

- **Bioconductor Guidelines**: https://contributions.bioconductor.org/
- **Package Guidelines**: https://contributions.bioconductor.org/package-guidelines.html
- **Submission Tracker**: https://github.com/Bioconductor/Contributions/issues
- **BiocCheck**: https://bioconductor.org/packages/BiocCheck/

## Support

- **Bioc-devel mailing list**: https://stat.ethz.ch/mailman/listinfo/bioc-devel
- **Slack**: https://bioc-community.herokuapp.com/
- **Support site**: https://support.bioconductor.org/

## Current Package Status

✅ **Ready for submission after completing Pre-Submission Checklist items**

**Statistics:**
- R CMD check: 0 errors, 0 warnings, 3 acceptable notes
- Test coverage: 25 tests, 100% passing
- Documentation: 34 help files, complete vignette
- Code quality: All functions documented, proper imports

**Next immediate steps:**
1. Update author information in DESCRIPTION
2. Create GitHub repository
3. Change version to 0.99.0
4. Run BiocCheck
5. Submit issue on Bioconductor/Contributions
