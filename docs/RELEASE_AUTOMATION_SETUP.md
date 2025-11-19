# DSTL Release Automation - Complete Setup

## 🎯 What You Get

**Before:** Manual process
1. Build package locally
2. Create GitHub release  
3. Upload files manually
4. Write release notes
5. Copy/paste checksums

**After:** Fully automated
```bash
git tag v1.0.0 && git push origin v1.0.0
```
**Done!** ✨

---

## 📦 What's Included

### Scripts (in `scripts/`)
- ✅ `quick_release.sh` - Interactive release wizard
- ✅ `create_release_package.sh` - Build release packages
- ✅ `test_release_package.sh` - Validate packages before release

### GitHub Actions (in `.github/workflows/`)
- ✅ `auto-release.yml` - Automatic release on tag push
- ✅ `release-package.yml` - Upload packages to existing releases

### Documentation (in `docs/`)
- ✅ `AUTOMATED_RELEASES.md` - Complete automation guide
- ✅ `RELEASE_PACKAGES.md` - Package creation details
- ✅ `RELEASE_QUICK_REFERENCE.md` - Quick reference card

### Examples (in `examples/`)
- ✅ `FetchContent-Release.cmake` - How to consume release packages
- ✅ `CPM-Example.cmake` - Alternative using CPM.cmake

---

## 🚀 Quick Start

### Option 1: Use the Interactive Helper
```bash
./scripts/quick_release.sh
```
Follow the prompts - it will:
- Check your git status
- Validate version format
- Build and test package locally (optional)
- Create and push the tag
- Show you where to monitor progress

### Option 2: Manual (Power Users)
```bash
# Tag and push
git tag v1.2.3
git push origin v1.2.3

# Monitor at https://github.com/jpsilverberg/DSTL/actions
```

---

## 🔄 How It Works

```
┌─────────────────┐
│  Push v1.2.3    │
│  git tag        │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ GitHub Actions  │
│  Triggered      │
└────────┬────────┘
         │
         ├──→ Checkout code + submodules
         │
         ├──→ Build release package
         │       • dstl-1.2.3-headers.tar.gz
         │       • dstl-1.2.3-headers.zip
         │
         ├──→ Generate SHA256 checksums
         │
         ├──→ Create release notes
         │       • Usage examples
         │       • Checksums
         │       • Changelog
         │
         └──→ Create GitHub Release
                 • Upload all files
                 • Publish automatically
                 
┌─────────────────┐
│  Users can now  │
│  download!      │
└─────────────────┘
```

---

## 📋 Pre-Release Checklist

Before running `./scripts/quick_release.sh` or pushing a tag:

- [ ] All changes committed
- [ ] Tests pass locally
- [ ] No uncommitted changes (`git status`)
- [ ] On correct branch (usually `main`)
- [ ] Submodules up to date
- [ ] Version number chosen (semantic versioning: MAJOR.MINOR.PATCH)

---

## 🔧 First-Time Setup

### 1. Enable GitHub Actions Permissions

Go to your repository settings:
```
Settings → Actions → General → Workflow permissions
```

Select: **"Read and write permissions"**

Click: **"Save"**

This allows GitHub Actions to create releases and upload files.

### 2. Test Locally First (Optional but Recommended)

```bash
# Build a test package
./scripts/create_release_package.sh 1.0.0-test

# Verify it works
./scripts/test_release_package.sh 1.0.0-test

# Clean up
rm -rf build-release install-release dstl-*.tar.gz dstl-*.zip
```

### 3. Create Your First Release

```bash
./scripts/quick_release.sh
# Enter version: 1.0.0
# Follow prompts
```

### 4. Verify

- Go to https://github.com/jpsilverberg/DSTL/actions
- Watch the workflow run (takes ~2-3 minutes)
- Check https://github.com/jpsilverberg/DSTL/releases

---

## 📊 What Gets Released

### Included in Package:
✅ All headers from `include/`, `Bits/include/`  
✅ Component headers (Mutex, Logger, DataPools, etc.)  
✅ CMake configuration files  
✅ Minimal source files for compiled components  

### Excluded from Package:
❌ Build directories  
❌ Test suites  
❌ Documentation build files  
❌ IDE configurations  
❌ Development tools  
❌ Coverage reports  
❌ Git history  

**Result:** Package is typically **10-50x smaller** than full git clone!

---

## 🎓 Examples for Users

Once you've created a release, users can consume it like this:

### Simple Example
```cmake
cmake_minimum_required(VERSION 3.14)
project(MyApp)

include(FetchContent)
FetchContent_Declare(dstl
  URL https://github.com/jpsilverberg/DSTL/releases/download/v1.2.3/dstl-1.2.3-headers.tar.gz
)
FetchContent_MakeAvailable(dstl)
find_package(DSTL CONFIG REQUIRED COMPONENTS dmutex dpools)

add_executable(myapp main.cpp)
target_link_libraries(myapp PRIVATE dstl::dmutex dstl::dpools)
```

### With Integrity Check
```cmake
FetchContent_Declare(dstl
  URL https://github.com/jpsilverberg/DSTL/releases/download/v1.2.3/dstl-1.2.3-headers.tar.gz
  URL_HASH SHA256=abc123...  # From checksums.txt in release
)
```

---

## 🐛 Common Issues

### "Resource not accessible by integration"

**Cause:** GitHub Actions doesn't have permission to create releases

**Fix:** Settings → Actions → General → Workflow permissions → "Read and write"

### Tag already exists

**Fix:**
```bash
git tag -d v1.2.3
git push origin :refs/tags/v1.2.3
# Now create it again
```

### Package too large

**Check what's included:**
```bash
tar -tzf dstl-1.2.3-headers.tar.gz | less
```

**Fix:** Update exclusion patterns in `scripts/create_release_package.sh`

### Workflow doesn't trigger

**Check:**
- Tag format must be `v*.*.*` (e.g., `v1.2.3`)
- Workflow file must be on the default branch
- Actions must be enabled in repository settings

---

## 📚 Learn More

| Document | Description |
|----------|-------------|
| [AUTOMATED_RELEASES.md](AUTOMATED_RELEASES.md) | Complete automation guide |
| [RELEASE_PACKAGES.md](RELEASE_PACKAGES.md) | Package creation details |
| [RELEASE_QUICK_REFERENCE.md](RELEASE_QUICK_REFERENCE.md) | Quick reference card |
| [README.md](../README.md) | Main project documentation |

---

## 🎉 You're All Set!

Your DSTL repository is now configured for fully automated releases!

**Next steps:**
1. Make sure GitHub Actions permissions are set
2. Run `./scripts/quick_release.sh` when ready
3. Let automation handle the rest

**Questions?** Check the documentation links above.

**Happy releasing! 🚀**
