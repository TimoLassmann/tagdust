# CMake Build System for TagDust2

This document describes the new CMake build system for TagDust2, which provides a modern alternative to the traditional GNU Autotools build system.

## Quick Start

### Building with CMake

```bash
# Create a build directory
mkdir build
cd build

# Configure the build
cmake ..

# Build the project  
make -j$(nproc)

# Run tests
ctest

# Install (optional)
make install
```

### Building with different configurations

```bash
# Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build (default)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Enable debugging with verbosity level 2
cmake -DENABLE_DEBUGGING=ON -DDEBUG_LEVEL=2 ..

# Enable Valgrind tests
cmake -DENABLE_VALGRIND_TESTS=ON ..

# Disable tests
cmake -DBUILD_TESTS=OFF ..
```

## Configuration Options

### Build Types
- **Release** (default): Optimized build with `-O2 -funroll-loops -Wall -std=gnu99`
- **Debug**: Debug build with `-ggdb -Wall -m64 -std=gnu99`
- **RelWithDebInfo**: Release with debug info
- **MinSizeRel**: Minimal size release

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `ENABLE_DEBUGGING` | OFF | Enable debugging code |
| `DEBUG_LEVEL` | 0 | Debug verbosity level (0-3) |
| `ENABLE_VALGRIND_TESTS` | OFF | Run tests under Valgrind |
| `BUILD_TESTS` | ON | Build test programs |

### Examples

```bash
# Full debug build with maximum verbosity
cmake -DENABLE_DEBUGGING=ON -DDEBUG_LEVEL=3 -DCMAKE_BUILD_TYPE=Debug ..

# Release build without tests
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF ..

# Valgrind testing (requires valgrind installed)
cmake -DENABLE_VALGRIND_TESTS=ON ..
```

## Built Executables

The CMake build produces the following executables in the `src/` directory:

### Main Programs
- **tagdust**: Main demultiplexing tool
- **simreads**: Read simulation utility  
- **evalres**: Results evaluation tool
- **merge**: Sequence merging utility
- **rename_qiime**: QIIME format renaming tool

### Test Programs (when `BUILD_TESTS=ON`)
- **tagdustiotest**: I/O unit tests
- **tagdust_rtest**: Reproducible tagdust with fixed seed
- **simreads_rtest**: Reproducible simreads with fixed seed  
- **evalres_rtest**: Reproducible evalres with fixed seed

## Testing

The CMake build includes comprehensive testing support:

```bash
# Run all tests
ctest

# Run tests with output
ctest --output-on-failure

# Run only integration tests
ctest -L integration

# Run a specific test
ctest -R sanity_test

# Run tests in parallel
ctest -j$(nproc)
```

### Test Categories
- **Unit tests**: `tagdustiotest` - Basic I/O functionality
- **Integration tests**: Shell script tests from `dev/` directory
  - `sanity_test`: Basic functionality tests
  - `bar_read_test`: Barcode reading tests
  - `casava_test`: CASAVA format tests

## Installation

```bash
# Install to system directories (default: /usr/local)
make install

# Install to custom location
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make install

# Install only binaries
make install/fast
```

### Install Components
- **Binaries**: Main executables in `${CMAKE_INSTALL_BINDIR}`
- **Headers**: Development headers in `${CMAKE_INSTALL_INCLUDEDIR}/tagdust2` (optional)
- **Documentation**: PDF and text docs in `${CMAKE_INSTALL_DOCDIR}`
- **Test data**: Test files and scripts in `${CMAKE_INSTALL_DATAROOTDIR}/tagdust2`

## Package Generation

```bash
# Create source packages
make package_source

# Available formats: TGZ, ZIP
```

## Dependencies

### Required
- CMake 3.18 or later
- C compiler with C11 support
- pthread library
- math library (libm)

### Optional
- **Pandoc**: For documentation generation
- **Valgrind**: For memory testing (if `ENABLE_VALGRIND_TESTS=ON`)

## Architecture Support

- **32-bit systems**: Supported
- **64-bit systems**: Supported with benchmark reproduction capability
- **SSE support**: Automatically detected (xmmintrin.h)

## Comparison with Autotools

### Advantages of CMake
- **Faster configuration**: No need for `./autogen.sh` step
- **Better IDE support**: Generates native IDE project files
- **Parallel builds**: Better parallel build support
- **Cross-platform**: Works on Windows, macOS, Linux
- **Modern syntax**: Cleaner, more readable build files
- **Better dependency handling**: Improved library detection

### Equivalent Commands

| Autotools | CMake |
|-----------|-------|
| `./autogen.sh && ./configure` | `cmake ..` |
| `make` | `make` or `cmake --build .` |
| `make check` | `ctest` |
| `make install` | `make install` |
| `./configure --enable-debugging` | `cmake -DENABLE_DEBUGGING=ON ..` |
| `./configure --enable-valgrind-tests` | `cmake -DENABLE_VALGRIND_TESTS=ON ..` |

## Troubleshooting

### Common Issues

1. **CMake too old**: Ensure CMake 3.18 or later
   ```bash
   cmake --version
   ```

2. **Missing dependencies**: Install required libraries
   ```bash
   # Ubuntu/Debian
   sudo apt-get install build-essential cmake libpthread-stubs0-dev
   
   # macOS
   brew install cmake
   ```

3. **Test failures**: Ensure executables are built
   ```bash
   ls -la src/tagdust*
   ```

4. **Permission errors**: Use proper build directory
   ```bash
   mkdir build && cd build  # Don't build in source directory
   ```

### Debug Build Issues

```bash
# Clean rebuild
rm -rf build
mkdir build && cd build
cmake .. && make
```

## Migration from Autotools

To migrate an existing autotools build:

1. Clean existing build artifacts:
   ```bash
   make distclean  # or git clean -fdx
   ```

2. Create CMake build:
   ```bash
   mkdir build && cd build
   cmake ..
   make
   ```

3. Verify equivalence:
   ```bash
   ctest  # Should pass all tests
   ```

## Future Enhancements

Potential improvements for the CMake build system:

- **SIMD optimization**: Automatic vectorization flags
- **Static analysis**: Integration with clang-tidy, cppcheck
- **Code coverage**: Coverage reporting with gcov/lcov
- **Continuous integration**: GitHub Actions integration
- **Documentation**: Automatic API documentation generation

## Support

For issues with the CMake build system, please:

1. Check this documentation
2. Verify CMake version and dependencies
3. Try a clean rebuild
4. Report issues with full error messages and system information

The CMake build system maintains full compatibility with the original autotools functionality while providing modern build features and improved maintainability.