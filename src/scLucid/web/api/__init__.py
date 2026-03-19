"""
FastAPI backend for scLucid web application.

Provides REST API for:
- Data upload and management
- QC analysis
- Preprocessing
- Clustering and annotation
- Visualization
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from scLucid.web.api.routes import qc, preprocess, analysis, data

app = FastAPI(
    title="scLucid API",
    description="Single-cell RNA-seq analysis API",
    version="0.1.0",
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(qc.router, prefix="/api/qc", tags=["QC"])
app.include_router(preprocess.router, prefix="/api/preprocess", tags=["Preprocessing"])
app.include_router(analysis.router, prefix="/api/analysis", tags=["Analysis"])
app.include_router(data.router, prefix="/api/data", tags=["Data"])


@app.get("/")
async def root():
    """Root endpoint."""
    return {
        "message": "scLucid API",
        "version": "0.1.0",
        "docs": "/docs",
    }


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy"}
