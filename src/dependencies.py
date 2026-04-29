"""
Módulo para verificação e instalação automática de dependências externas.
Garante que BLAST+, BEDTools e SSEARCH estejam disponíveis.
"""

import subprocess
import logging
import os
import multiprocessing
from pathlib import Path
from typing import Tuple, Optional

logger = logging.getLogger(__name__)

# Diretório base do projeto
PROJECT_DIR = Path(__file__).parent.parent

# Diretórios das ferramentas
BLAST_DIR = PROJECT_DIR / "BLAST" / "ncbi-blast-2.17.0+"
BEDTOOLS_DIR = PROJECT_DIR / "BEDTOOLS" / "bedtools2"
FASTA36_DIR = PROJECT_DIR / "FASTA36" / "fasta-36.3.8i"

# URLs de download
BLAST_URL = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz"
BEDTOOLS_URL = "https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz"


def get_num_threads() -> int:
    """
    Retorna o número de threads a ser usado (total de CPUs - 1).
    Mínimo de 1 thread.
    
    Returns:
        Número de threads
    """
    num_cpus = multiprocessing.cpu_count()
    num_threads = max(1, num_cpus - 1)
    logger.debug(f"CPUs disponíveis: {num_cpus}, usando {num_threads} threads")
    return num_threads


def check_tool_in_path(tool_name: str) -> bool:
    """
    Verifica se uma ferramenta está disponível no PATH do sistema.
    
    Args:
        tool_name: Nome da ferramenta
    
    Returns:
        True se encontrada, False caso contrário
    """
    try:
        result = subprocess.run(
            ["which", tool_name],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except Exception:
        return False


def install_blast() -> Path:
    """
    Instala BLAST+ localmente se não estiver disponível.
    
    Returns:
        Caminho para o diretório bin do BLAST
    
    Raises:
        RuntimeError: Se falhar a instalação
    """
    logger.info("Instalando BLAST+ localmente...")
    
    blast_parent = PROJECT_DIR / "BLAST"
    blast_parent.mkdir(parents=True, exist_ok=True)
    
    if BLAST_DIR.exists():
        logger.info(f"BLAST+ já instalado em: {BLAST_DIR}")
        return BLAST_DIR / "bin"
    
    try:
        # Download
        logger.info(f"Baixando BLAST+ de: {BLAST_URL}")
        tar_file = blast_parent / "blast.tar.gz"
        
        subprocess.run(
            ["wget", "-q", "-O", str(tar_file), BLAST_URL],
            check=True
        )
        
        # Extração
        logger.info("Extraindo BLAST+...")
        subprocess.run(
            ["tar", "-xzf", str(tar_file), "-C", str(blast_parent)],
            check=True
        )
        
        # Limpeza
        tar_file.unlink()
        
        logger.info(f"BLAST+ instalado com sucesso em: {BLAST_DIR}")
        return BLAST_DIR / "bin"
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Erro ao instalar BLAST+: {e}")
    except Exception as e:
        raise RuntimeError(f"Erro inesperado ao instalar BLAST+: {e}")


def install_bedtools() -> Path:
    """
    Instala BEDTools localmente se não estiver disponível.
    
    Returns:
        Caminho para o diretório bin do BEDTools
    
    Raises:
        RuntimeError: Se falhar a instalação
    """
    logger.info("Instalando BEDTools localmente...")
    
    bedtools_parent = PROJECT_DIR / "BEDTOOLS"
    bedtools_parent.mkdir(parents=True, exist_ok=True)
    
    if BEDTOOLS_DIR.exists() and (BEDTOOLS_DIR / "bin" / "bedtools").exists():
        logger.info(f"BEDTools já instalado em: {BEDTOOLS_DIR}")
        return BEDTOOLS_DIR / "bin"
    
    try:
        # Download
        logger.info(f"Baixando BEDTools de: {BEDTOOLS_URL}")
        tar_file = bedtools_parent / "bedtools.tar.gz"
        
        subprocess.run(
            ["wget", "-q", "-O", str(tar_file), BEDTOOLS_URL],
            check=True
        )
        
        # Extração
        logger.info("Extraindo BEDTools...")
        subprocess.run(
            ["tar", "-xzf", str(tar_file), "-C", str(bedtools_parent)],
            check=True
        )
        
        # Compilação
        logger.info("Compilando BEDTools (isso pode levar alguns minutos)...")
        num_threads = get_num_threads()
        
        subprocess.run(
            ["make", f"-j{num_threads}"],
            cwd=str(BEDTOOLS_DIR),
            check=True,
            capture_output=True
        )
        
        # Limpeza
        tar_file.unlink()
        
        logger.info(f"BEDTools instalado com sucesso em: {BEDTOOLS_DIR}")
        return BEDTOOLS_DIR / "bin"
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Erro ao instalar BEDTools: {e}")
    except Exception as e:
        raise RuntimeError(f"Erro inesperado ao instalar BEDTools: {e}")


def get_blast_paths() -> Tuple[str, str]:
    """
    Retorna os caminhos para makeblastdb e tblastn.
    FORÇA uso da instalação local para evitar problemas com versões do sistema.
    
    Returns:
        Tupla (makeblastdb_path, tblastn_path)
    """
    # FORÇA instalação local (ignora sistema)
    if BLAST_DIR.exists():
        makeblastdb_path = str(BLAST_DIR / "bin" / "makeblastdb")
        tblastn_path = str(BLAST_DIR / "bin" / "tblastn")
        logger.debug("Using LOCAL BLAST+ (forced)")
        return makeblastdb_path, tblastn_path
    
    # Instala localmente
    logger.info("Installing BLAST+ locally (forced local installation)...")
    blast_bin = install_blast()
    return str(blast_bin / "makeblastdb"), str(blast_bin / "tblastn")


def get_bedtools_path() -> str:
    """
    Return path to bedtools.
    FORÇA uso da instalação local para evitar problemas com versões do sistema.
    
    Returns:
        Path to bedtools executable
    """
    # FORÇA instalação local (ignora sistema)
    if BEDTOOLS_DIR.exists() and (BEDTOOLS_DIR / "bin" / "bedtools").exists():
        logger.debug("Using LOCAL BEDTools (forced)")
        return str(BEDTOOLS_DIR / "bin" / "bedtools")
    
    # Install locally
    logger.info("Installing BEDTools locally (forced local installation)...")
    bedtools_bin = install_bedtools()
    return str(bedtools_bin / "bedtools")


def install_fasta36() -> Path:
    """
    Instala FASTA36 suite localmente se não estiver disponível.
    
    Returns:
        Caminho para o diretório bin do FASTA36
    """
    logger.debug("Installing FASTA36 locally...")
    
    fasta36_url = "https://fasta.bioch.virginia.edu/wrpearson/fasta/CURRENT/fasta36-linux64.tar.gz"
    
    try:
        # Cria diretório
        FASTA36_DIR.parent.mkdir(parents=True, exist_ok=True)
        
        # Download do arquivo
        tar_file = FASTA36_DIR.parent / "fasta36-linux64.tar.gz"
        logger.info(f"Baixando FASTA36 de {fasta36_url}...")
        
        subprocess.run(
            ["wget", "-O", str(tar_file), fasta36_url],
            check=True,
            capture_output=True
        )
        
        # Extrai o arquivo
        logger.info("Extraindo FASTA36...")
        subprocess.run(
            ["tar", "-xzf", str(tar_file), "-C", str(FASTA36_DIR.parent)],
            check=True,
            capture_output=True
        )
        
        # Remove o tar.gz
        tar_file.unlink()
        
        # Torna executáveis
        bin_dir = FASTA36_DIR / "bin"
        if bin_dir.exists():
            for exe in bin_dir.glob("*"):
                if exe.is_file():
                    exe.chmod(0o755)
        
        logger.info(f"✓ FASTA36 instalado com sucesso em: {FASTA36_DIR}")
        return bin_dir
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Erro ao instalar FASTA36: {e}")
    except Exception as e:
        raise RuntimeError(f"Erro inesperado ao instalar FASTA36: {e}")


def get_ssearch_path() -> str:
    """
    Return path to ssearch36 (required).
    FORÇA uso da instalação local para evitar problemas com versões do sistema.
    
    Returns:
        Path to ssearch36
    
    Raises:
        RuntimeError: If ssearch36 not found
    """
    # FORÇA instalação local (ignora sistema)
    if FASTA36_DIR.exists() and (FASTA36_DIR / "bin" / "ssearch36").exists():
        ssearch_path = str(FASTA36_DIR / "bin" / "ssearch36")
        logger.debug(f"Using LOCAL SSEARCH36 (forced): {ssearch_path}")
        return ssearch_path
    
    # SSEARCH36 not found - install automatically
    logger.info("Installing FASTA36 locally (forced local installation)...")
    install_fasta36()
    
    if FASTA36_DIR.exists() and (FASTA36_DIR / "bin" / "ssearch36").exists():
        ssearch_path = str(FASTA36_DIR / "bin" / "ssearch36")
        logger.info(f"✓ SSEARCH36 instalado com sucesso: {ssearch_path}")
        return ssearch_path
    
    raise RuntimeError(
        "SSEARCH36 não encontrado e falha na instalação automática.\n"
        f"Esperado em: {FASTA36_DIR / 'bin' / 'ssearch36'}"
    )


def verify_and_install_dependencies() -> dict:
    """
    Check and install all required dependencies.
    
    Returns:
        Dictionary with tool paths
    """
    dependencies = {}
    
    # BLAST+ (silent installation)
    makeblastdb_path, tblastn_path = get_blast_paths()
    dependencies['makeblastdb'] = makeblastdb_path
    dependencies['tblastn'] = tblastn_path
    
    # BEDTools (silent installation)
    bedtools_path = get_bedtools_path()
    dependencies['bedtools'] = bedtools_path
    
    # SSEARCH36 (silent installation)
    ssearch_path = get_ssearch_path()
    dependencies['ssearch36'] = ssearch_path
    
    # Thread count
    num_threads = get_num_threads()
    dependencies['num_threads'] = num_threads
    
    logger.info(f"Dependencies verified ({num_threads} threads available)")
    
    return dependencies


if __name__ == '__main__':
    verify_and_install_dependencies()