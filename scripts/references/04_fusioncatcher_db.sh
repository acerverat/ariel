#!/bin/bash
set -e

# Base de datos de FusionCatcher
# El script de descarga se extrae del contenedor en tiempo de ejecucion para
# que siempre este actualizado. Se ejecuta en el HOST (no dentro de Docker)
# para que las descargas parciales persistan entre reinicios. Un wrapper de
# wget inyectado via PATH agrega reanudacion (-c) y reintentos sin tocar el
# script original de FusionCatcher.
#
# En tiempo de ejecucion del pipeline, fusioncatcher_db/ se monta en
# /opt/fusioncatcher/data dentro del contenedor.

rutaReferencias="$1"
outdir="$rutaReferencias/GRCh38_no_alt"
dbdir="$outdir/fusioncatcher_db"
sentinel="$dbdir/current/version.txt"

if [[ -f "$sentinel" ]]; then
    echo "[fusioncatcher_db] Base de datos ya presente, omitiendo descarga."
    exit 0
fi

mkdir -p "$dbdir"

# Extraer el script de descarga directamente del contenedor para que
# cualquier actualizacion de FusionCatcher se refleje automaticamente.
DOWNLOAD_SCRIPT="$dbdir/download-human-db.sh"
echo "[fusioncatcher_db] Extrayendo script de descarga del contenedor..."
docker run --rm acerverat/ariel-env:latest \
    cat /opt/fusioncatcher/data/download-human-db.sh > "$DOWNLOAD_SCRIPT"
chmod +x "$DOWNLOAD_SCRIPT"

# Wrapper de wget que agrega reanudacion y rotacion de mirrors de SourceForge
# sin modificar el script original. Se coloca primero en PATH para que el
# script lo use en lugar del wget del sistema.
#
# Cuando detecta una URL de SourceForge (original o de mirror), la convierte
# en una URL directa de mirror y rota por la lista hasta encontrar uno que
# funcione. Cada mirror se intenta con -c para reanudar descargas parciales.
WGET_WRAPPER_DIR="$(mktemp -d)"
cat > "$WGET_WRAPPER_DIR/wget" << 'EOF'
#!/bin/bash

SF_MIRRORS=(
    "netix.dl.sourceforge.net"
    "netcologne.dl.sourceforge.net"
    "jaist.dl.sourceforge.net"
    "cfhcable.dl.sourceforge.net"
    "cytranet.dl.sourceforge.net"
    "astuteinternet.dl.sourceforge.net"
    "vorboss.dl.sourceforge.net"
    "excellmedia.dl.sourceforge.net"
    "phoenixnap.dl.sourceforge.net"
)

# Busca si algun argumento es una URL de SourceForge
IS_SF=0
for arg in "$@"; do
    if [[ "$arg" == *"sourceforge.net"* ]]; then
        IS_SF=1
        break
    fi
done

if [[ "$IS_SF" -eq 0 ]]; then
    exec /usr/bin/wget --tries=50 --waitretry=30 --retry-connrefused \
                       --timeout=60 -c "$@"
fi

# Extrae la ruta del archivo de la URL (funciona con URLs originales de SF
# y con URLs directas de mirrors).
FILE_PATH=""
for arg in "$@"; do
    if [[ "$arg" == *"sourceforge.net/projects/fusioncatcher/files/data/"* ]]; then
        FILE_PATH="${arg##*/projects/fusioncatcher/files/data/}"
        FILE_PATH="${FILE_PATH%%\?*}"
        break
    elif [[ "$arg" == *".dl.sourceforge.net/project/fusioncatcher/data/"* ]]; then
        FILE_PATH="${arg##*/project/fusioncatcher/data/}"
        FILE_PATH="${FILE_PATH%%\?*}"
        break
    fi
done

if [[ -z "$FILE_PATH" ]]; then
    exec /usr/bin/wget --tries=50 --waitretry=30 --retry-connrefused \
                       --timeout=60 -c "$@"
fi

# Reemplaza la URL en los args con la del mirror actual e intenta la descarga
try_mirror() {
    local mirror="$1"
    local mirror_url="https://${mirror}/project/fusioncatcher/data/${FILE_PATH}"
    local new_args=()
    for arg in "$@"; do
        if [[ "$arg" == *"sourceforge.net"* ]]; then
            new_args+=("$mirror_url")
        else
            new_args+=("$arg")
        fi
    done
    echo "[wget-wrapper] Mirror: $mirror" >&2
    /usr/bin/wget --tries=3 --waitretry=15 --retry-connrefused \
                  --timeout=60 -c "${new_args[@]}"
}

for mirror in "${SF_MIRRORS[@]}"; do
    if try_mirror "$mirror" "$@"; then
        exit 0
    fi
    echo "[wget-wrapper] Mirror $mirror fallo, probando el siguiente..." >&2
done

echo "[wget-wrapper] Todos los mirrors fallaron." >&2
exit 1
EOF
chmod +x "$WGET_WRAPPER_DIR/wget"

echo "[fusioncatcher_db] Descargando base de datos de FusionCatcher..."
cd "$dbdir"
PATH="$WGET_WRAPPER_DIR:$PATH" bash "$DOWNLOAD_SCRIPT"

rm -rf "$WGET_WRAPPER_DIR"
echo "[fusioncatcher_db] Listo."
