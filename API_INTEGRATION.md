# API Integration Guide for DESandbox

## Overview

DESandbox is designed to integrate seamlessly with web applications and REST APIs. This guide demonstrates how to use DESandbox as a backend for Django, Flask, or other web frameworks.

## Architecture Pattern

```
Frontend (React/Vue) 
    ↓ HTTP Request (JSON)
API Layer (Django/Flask)
    ↓ Rscript/rpy2
DESandbox R Package
    ↓ Analysis Results
API Layer (JSON Response)
    ↓
Frontend (Visualization)
```

## Integration Approaches

### 1. Command-Line Interface (Recommended for Django)

Create an R script that accepts JSON input and produces JSON output.

#### R Script: `desandbox_api.R`

```r
#!/usr/bin/env Rscript

library(DESandbox)
library(jsonlite)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript desandbox_api.R <input.json> <output.json>")
}

input_file <- args[1]
output_file <- args[2]

# Read input JSON
input_data <- fromJSON(input_file)

# Extract data
counts <- as.matrix(input_data$counts)
metadata <- as.data.frame(input_data$metadata)
params <- input_data$parameters

# Run analysis
tryCatch({
  # Create DESandbox object
  dso <- create_desandbox_object(
    counts = counts,
    metadata = metadata,
    condition_column = params$condition_column
  )
  
  # Run selected methods
  methods <- params$methods %||% c("DESeq2", "edgeR", "limma-voom")
  results <- run_all_methods(dso, methods = methods)
  
  # Filter results
  filtered_results <- lapply(results, function(res) {
    filter_results(
      res$results,
      padj_threshold = params$padj_threshold %||% 0.05,
      lfc_threshold = params$lfc_threshold %||% 1
    )
  })
  
  # Compare methods
  comparison <- compare_methods(
    results,
    padj_threshold = params$padj_threshold %||% 0.05,
    lfc_threshold = params$lfc_threshold %||% 1
  )
  
  # Prepare output
  output <- list(
    status = "success",
    results = filtered_results,
    comparison = list(
      summary_stats = comparison$summary_stats,
      overlap_matrix = comparison$comparison_matrix,
      consensus_genes = comparison$consensus_genes
    ),
    timestamp = Sys.time()
  )
  
}, error = function(e) {
  output <<- list(
    status = "error",
    message = e$message,
    timestamp = Sys.time()
  )
})

# Write output JSON
write(toJSON(output, auto_unbox = TRUE, pretty = TRUE), output_file)
```

#### Django View: `views.py`

```python
import json
import subprocess
import tempfile
import pandas as pd
from pathlib import Path
from django.http import JsonResponse
from django.views import View
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status


class DifferentialExpressionView(APIView):
    """
    Run differential expression analysis using DESandbox
    """
    
    def post(self, request):
        """
        Expected JSON format:
        {
            "counts": [[...], [...]],  // Gene x Sample matrix
            "metadata": {
                "sample": [...],
                "condition": [...]
            },
            "parameters": {
                "condition_column": "condition",
                "methods": ["DESeq2", "edgeR", "limma-voom"],
                "padj_threshold": 0.05,
                "lfc_threshold": 1
            }
        }
        """
        
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as input_file:
                json.dump(request.data, input_file)
                input_path = input_file.name
            
            output_path = tempfile.mktemp(suffix='.json')
            
            # Run R script
            result = subprocess.run(
                ['Rscript', 'desandbox_api.R', input_path, output_path],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            # Check for errors
            if result.returncode != 0:
                return Response(
                    {'error': 'R script failed', 'stderr': result.stderr},
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR
                )
            
            # Read output
            with open(output_path, 'r') as f:
                output_data = json.load(f)
            
            # Clean up temporary files
            Path(input_path).unlink()
            Path(output_path).unlink()
            
            if output_data.get('status') == 'error':
                return Response(
                    {'error': output_data.get('message')},
                    status=status.HTTP_400_BAD_REQUEST
                )
            
            return Response(output_data, status=status.HTTP_200_OK)
            
        except subprocess.TimeoutExpired:
            return Response(
                {'error': 'Analysis timeout'},
                status=status.HTTP_408_REQUEST_TIMEOUT
            )
        except Exception as e:
            return Response(
                {'error': str(e)},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )


class EnrichmentAnalysisView(APIView):
    """
    Run enrichment analysis on gene list
    """
    
    def post(self, request):
        """
        Expected JSON format:
        {
            "genes": ["GENE1", "GENE2", ...],
            "organism": "human",
            "analysis_type": "both",
            "ontology": "BP"
        }
        """
        
        # Similar pattern as above
        pass
```

### 2. Using rpy2 (Direct Python-R Integration)

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Enable pandas-R conversion
pandas2ri.activate()

# Import DESandbox
desandbox = importr('DESandbox')

class DESandboxService:
    """
    Service class for DESandbox integration
    """
    
    def __init__(self):
        self.desandbox = importr('DESandbox')
    
    def run_analysis(self, counts_df, metadata_df, condition_column='condition'):
        """
        Run differential expression analysis
        
        Args:
            counts_df: pandas DataFrame (genes x samples)
            metadata_df: pandas DataFrame with sample metadata
            condition_column: name of condition column
        
        Returns:
            dict: Analysis results
        """
        
        # Convert pandas to R
        r_counts = pandas2ri.py2rpy(counts_df)
        r_metadata = pandas2ri.py2rpy(metadata_df)
        
        # Create DESandbox object
        dso = self.desandbox.create_desandbox_object(
            counts=r_counts,
            metadata=r_metadata,
            condition_column=condition_column
        )
        
        # Run all methods
        results = self.desandbox.run_all_methods(dso)
        
        # Convert results back to Python
        py_results = {}
        for method in results.names:
            method_result = results.rx2(method)
            py_results[method] = {
                'results': pandas2ri.rpy2py(method_result.rx2('results')),
                'parameters': dict(method_result.rx2('parameters'))
            }
        
        return py_results
    
    def filter_results(self, results_df, padj_threshold=0.05, lfc_threshold=1):
        """
        Filter DE results
        """
        r_results = pandas2ri.py2rpy(results_df)
        filtered = self.desandbox.filter_results(
            r_results,
            padj_threshold=padj_threshold,
            lfc_threshold=lfc_threshold
        )
        return pandas2ri.rpy2py(filtered)


# Usage in Django view
service = DESandboxService()

class DifferentialExpressionView(APIView):
    def post(self, request):
        counts_df = pd.DataFrame(request.data['counts'])
        metadata_df = pd.DataFrame(request.data['metadata'])
        
        results = service.run_analysis(counts_df, metadata_df)
        
        # Convert to JSON-serializable format
        response_data = {
            method: {
                'results': res['results'].to_dict(orient='records'),
                'parameters': res['parameters']
            }
            for method, res in results.items()
        }
        
        return Response(response_data)
```

### 3. REST API with Plumber (R-based API)

Create an R-based REST API using the `plumber` package:

#### `api.R`

```r
library(plumber)
library(DESandbox)
library(jsonlite)

#* @apiTitle DESandbox REST API
#* @apiDescription Differential expression analysis API

#* Run differential expression analysis
#* @param counts:object Count matrix (JSON)
#* @param metadata:object Sample metadata (JSON)
#* @param condition_column:string Condition column name
#* @param methods:[string] Methods to run
#* @param padj_threshold:numeric Adjusted p-value threshold
#* @param lfc_threshold:numeric Log2 fold change threshold
#* @post /analyze
#* @serializer json
function(counts, metadata, condition_column = "condition", 
         methods = c("DESeq2", "edgeR", "limma-voom"),
         padj_threshold = 0.05, lfc_threshold = 1) {
  
  tryCatch({
    # Parse JSON input
    counts_mat <- as.matrix(counts)
    metadata_df <- as.data.frame(metadata)
    
    # Create DESandbox object
    dso <- create_desandbox_object(
      counts = counts_mat,
      metadata = metadata_df,
      condition_column = condition_column
    )
    
    # Run analysis
    results <- run_all_methods(dso, methods = methods)
    
    # Filter results
    filtered <- lapply(results, function(res) {
      filter_results(
        res$results,
        padj_threshold = padj_threshold,
        lfc_threshold = lfc_threshold
      )
    })
    
    # Compare methods
    comparison <- compare_methods(results, padj_threshold, lfc_threshold)
    
    list(
      status = "success",
      results = filtered,
      comparison = comparison$summary_stats
    )
    
  }, error = function(e) {
    list(
      status = "error",
      message = e$message
    )
  })
}

#* Run enrichment analysis
#* @param genes:[string] Gene list
#* @param organism:string Organism (human, mouse, rat)
#* @param analysis_type:string Analysis type (GO, KEGG, both)
#* @post /enrichment
#* @serializer json
function(genes, organism = "human", analysis_type = "both") {
  
  tryCatch({
    enrichment <- run_enrichment(
      gene_list = genes,
      organism = organism,
      analysis_type = analysis_type
    )
    
    list(
      status = "success",
      GO = if (!is.null(enrichment$GO)) enrichment$GO$summary else NULL,
      KEGG = if (!is.null(enrichment$KEGG)) enrichment$KEGG$summary else NULL
    )
    
  }, error = function(e) {
    list(
      status = "error",
      message = e$message
    )
  })
}

#* Health check
#* @get /health
function() {
  list(status = "ok", timestamp = Sys.time())
}
```

#### Start the API server:

```r
library(plumber)
pr <- plumb("api.R")
pr$run(host = "0.0.0.0", port = 8000)
```

#### Django integration:

```python
import requests

class DESandboxAPIClient:
    def __init__(self, base_url="http://localhost:8000"):
        self.base_url = base_url
    
    def run_analysis(self, counts, metadata, **kwargs):
        response = requests.post(
            f"{self.base_url}/analyze",
            json={
                "counts": counts,
                "metadata": metadata,
                **kwargs
            }
        )
        return response.json()
    
    def run_enrichment(self, genes, organism="human"):
        response = requests.post(
            f"{self.base_url}/enrichment",
            json={"genes": genes, "organism": organism}
        )
        return response.json()
```

## Data Format Specifications

### Input Format: Counts Matrix

```json
{
  "counts": {
    "gene_names": ["Gene1", "Gene2", "Gene3"],
    "sample_names": ["Sample1", "Sample2", "Sample3", "Sample4"],
    "values": [
      [100, 150, 200, 180],
      [50, 60, 55, 58],
      [300, 320, 310, 305]
    ]
  }
}
```

Or as a flat structure:

```json
{
  "counts": [
    {"gene": "Gene1", "Sample1": 100, "Sample2": 150, ...},
    {"gene": "Gene2", "Sample1": 50, "Sample2": 60, ...}
  ]
}
```

### Input Format: Metadata

```json
{
  "metadata": {
    "sample": ["Sample1", "Sample2", "Sample3", "Sample4"],
    "condition": ["control", "control", "treatment", "treatment"],
    "batch": ["A", "B", "A", "B"]
  }
}
```

### Output Format: Results

```json
{
  "status": "success",
  "results": {
    "DESeq2": [
      {
        "gene": "Gene1",
        "baseMean": 157.5,
        "log2FoldChange": 2.34,
        "pvalue": 0.001,
        "padj": 0.01
      }
    ],
    "edgeR": [...],
    "limma-voom": [...]
  },
  "comparison": {
    "summary_stats": {
      "method": ["DESeq2", "edgeR", "limma-voom"],
      "n_significant": [120, 115, 118],
      "n_upregulated": [65, 60, 62],
      "n_downregulated": [55, 55, 56]
    },
    "consensus_genes": ["Gene1", "Gene5", "Gene10"]
  },
  "timestamp": "2026-01-02T10:30:00Z"
}
```

## Deployment Considerations

### Docker Container

```dockerfile
FROM rocker/r-ver:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install R packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'clusterProfiler'))"
RUN R -e "install.packages(c('plumber', 'jsonlite', 'digest'))"

# Install DESandbox
COPY . /DESandbox
RUN R -e "install.packages('/DESandbox', repos = NULL, type = 'source')"

# Expose port
EXPOSE 8000

# Start API
CMD ["R", "-e", "pr <- plumber::plumb('/app/api.R'); pr$run(host='0.0.0.0', port=8000)"]
```

### Performance Optimization

1. **Caching**: Cache results by input hash
2. **Queue system**: Use Celery for long-running analyses
3. **Parallel processing**: Enable BiocParallel in R
4. **Resource limits**: Set memory/CPU limits per request

### Security Considerations

1. **Input validation**: Validate count matrix dimensions
2. **Rate limiting**: Prevent abuse
3. **Timeout**: Set maximum execution time
4. **Sandboxing**: Run R in isolated environment
5. **Authentication**: Require API keys for production

## Example: Complete Django Integration

See `examples/django_integration/` directory for a complete working example with:
- Django REST Framework views
- Celery task queue
- Redis caching
- React frontend
- Docker Compose setup

## Support

For integration questions, open an issue on GitHub or contact the maintainers.
