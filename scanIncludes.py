#!/usr/bin/env python3
"""
Include Scanner

This script scans all source files (CPP and header files) in the codebase and creates a map showing
which directories include headers from other directories. It then detects circular dependencies
between directories and reports which specific header files create each dependency.

The map format is: [directoryA, directoryB] -> [headerFileName]
where directoryA has a source file that includes a header from directoryB.

Usage:
    python scanIncludes.py [solution_dir]

If no solution_dir is provided, it defaults to the current directory + /src.
"""

import os
import sys
import re
from collections import defaultdict
from typing import Dict, List, Set, Tuple

class IncludeScanner:
    def __init__(self, solutionDir: str):
        self.solutionDir = solutionDir
        # Define all file extensions that should be scanned for #include statements
        self.sourceExtensions = {'.cpp', '.cc', '.cxx', '.h', '.hpp', '.hxx', '.hh'}
        # Map: (sourceDirectory, targetDirectory) -> set of header filenames
        self.includeMap: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        
    def isSourceFile(self, filename: str) -> bool:
        """Check if a file is a source file based on its extension."""
        _, ext = os.path.splitext(filename)
        return ext.lower() in self.sourceExtensions
    
    def extractIncludesFromSource(self, sourceFilePath: str, sourceDir: str) -> Set[Tuple[str, str]]:
        """
        Extract all #include statements from a source file and return (headerName, headerDir) tuples.
        
        Only processes quoted includes (#include "file.h"), not angle bracket includes (#include <file.h>).
        Handles both relative paths and same-directory includes.
        """
        includes = set()
        
        try:
            with open(sourceFilePath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Find all #include statements with quoted includes only
                # This regex matches only #include "file" (not #include <file>)
                includePattern = r'#include\s*"([^"]+)"'
                matches = re.findall(includePattern, content)
                
                for includePath in matches:
                    # Determine header location based on path
                    if '/' in includePath or '\\' in includePath:
                        # Path contains slashes - relative to solution directory
                        headerDir = os.path.dirname(includePath)
                        headerName = os.path.basename(includePath)
                        
                        # Verify the header exists at this location using full path
                        # includePath is relative to the solution directory (src)
                        headerFullPath = os.path.abspath(os.path.join(self.solutionDir, includePath))
                        if os.path.exists(headerFullPath):
                            includes.add((headerName, headerDir))
                        else:
                            print(f"ERROR: Header not found at {headerFullPath}")
                            print(f"  Referenced in: {sourceFilePath}")
                            print(f"  Include path: {includePath}")
                            print(f"  Solution dir: {self.solutionDir}")
                            sys.exit(1)
                    else:
                        # No slashes - header is in same directory as source file
                        headerName = includePath
                        headerDir = sourceDir
                        
                        # Verify the header exists in the same directory using full path
                        headerFullPath = os.path.abspath(os.path.join(self.solutionDir, sourceDir, headerName))
                        if os.path.exists(headerFullPath):
                            # Don't add to includes since it's in the same directory
                            pass
                        else:
                            print(f"ERROR: Header not found at {headerFullPath}")
                            print(f"  Referenced in: {sourceFilePath}")
                            print(f"  Include path: {includePath}")
                            print(f"  Source dir: {sourceDir}")
                            sys.exit(1)
                    
        except Exception as e:
            print(f"Warning: Could not read {sourceFilePath}: {e}")
        
        return includes
    
    def scanSourceFiles(self):
        """
        Scan all source files (CPP and headers) and build the include map.
        
        Walks through the solution directory, finds all source files, extracts their
        #include statements, and builds a map showing which directories depend on
        headers from other directories.
        """
        print(f"Scanning source files in: {os.path.abspath(self.solutionDir)}")
        
        sourceFilesFound = 0
        
        for root, dirs, files in os.walk(self.solutionDir):
            # Skip common directories that shouldn't contain source files
            dirs[:] = [d for d in dirs if not d.startswith('.') and d not in {'build', 'bin', 'obj', 'Debug', 'Release', 'x64', 'x86', 'external'}]
            
            for file in files:
                if self.isSourceFile(file):
                    sourceFilesFound += 1
                    sourceFilePath = os.path.join(root, file)
                    sourceDir = os.path.relpath(root, self.solutionDir)
                    
                    # Extract includes from this source file
                    includes = self.extractIncludesFromSource(sourceFilePath, sourceDir)
                    
                    # For each included header, create the map entry
                    for headerName, headerDir in includes:
                        if sourceDir != headerDir:  # Only include if different directories
                            self.includeMap[(sourceDir, headerDir)].add(headerName)
        
        print(f"Scanned {sourceFilesFound} source files")
    
    def checkCircularDependencies(self):
        """
        Check for circular dependencies in the include map using depth-first search.
        
        Uses DFS to detect cycles in the dependency graph. When a cycle is found,
        stores the cycle path for later reporting.
        
        Returns:
            bool: True if a circular dependency is found, False otherwise
        """
        self.cyclePath = []
        
        def hasCycle(node, path, visited):
            """
            Recursive DFS function to detect cycles.
            
            Args:
                node: Current directory being visited
                path: Current path of directories visited
                visited: Set of directories already processed from this starting point
            
            Returns:
                bool: True if a cycle is found, False otherwise
            """
            if node in path:
                # Found a cycle - find the cycle path
                cycleStart = path.index(node)
                cycle = path[cycleStart:] + [node]
                self.cyclePath.extend(cycle)
                return True
            
            if node in visited:
                return False  # Already processed from this starting point
            
            visited.add(node)
            
            # Check all dependencies of this node
            for (sourceDir, targetDir), headers in self.includeMap.items():
                if sourceDir == node:
                    newPath = path.copy();
                    newPath.append(node);
                    if hasCycle(targetDir, newPath, visited):
                        return True
            
            return False
        
        # Check each directory for cycles
        allDirs = set()
        for (sourceDir, targetDir) in self.includeMap.keys():
            allDirs.add(sourceDir)
            allDirs.add(targetDir)
        
        for directory in allDirs:
            visited = set()  # Fresh visited set for each starting point
            if hasCycle(directory, [], visited):
                return True  # Stop at first cycle found
        
        return False

def main():
    """
    Main function that orchestrates the include scanning and circular dependency detection.
    
    Scans all source files, builds the dependency map, and reports any circular dependencies
    with details about which header files create each dependency in the cycle.
    """
    # Get solution directory from command line or use the directory where this script is located
    if len(sys.argv) > 1:
        solutionDir = sys.argv[1]
    else:
        # Use the directory where this script is located, plus /src
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        solutionDir = os.path.join(scriptDir, "src")
    
    if not os.path.exists(solutionDir):
        print(f"Error: Directory '{solutionDir}' does not exist.")
        sys.exit(1)
    
    scanner = IncludeScanner(solutionDir)
    
    # Scan source files and build the include map
    scanner.scanSourceFiles()
    
    # Check for circular dependencies
    if scanner.checkCircularDependencies():
        print("=" * 60)
        print("Circular dependency detected! Please restructure the code to break the circular dependency.")
        print("Cycle path:")
        
        # Iterate through the cycle path and show headers for each dependency
        for i in range(len(scanner.cyclePath) - 1):
            sourceDir = scanner.cyclePath[i]
            targetDir = scanner.cyclePath[i + 1]
            
            # Look up the headers that create this dependency
            headers = scanner.includeMap.get((sourceDir, targetDir), set())
            headerList = ", ".join(sorted(headers))
            
            print(f"  {sourceDir} -> {targetDir} (via: {headerList})")
        
        print("=" * 60)
        sys.exit(1)
    else:
        print(f"No circular dependencies detected.")

if __name__ == "__main__":
    main() 