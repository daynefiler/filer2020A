#!/bin/sh
TOOLDIR=$(pwd)/tools
BINDIR=${TOOLDIR}/bin
ERDS=7f47b69943e254e5ff5a21c7ab8915e9a4da568a

mkdir -p ${BINDIR}

wget "https://github.com/igm-team/ERDS/archive/${ERDS}.zip" -O "${TOOLDIR}/ERDS-${ERDS}.zip"
unzip ${TOOLDIR}/ERDS-${ERDS}.zip -d ${TOOLDIR}

pushd ${TOOLDIR}/ERDS-${ERDS}/erds_tcag/src/phmm; make; popd
pushd ${TOOLDIR}/ERDS-${ERDS}/erds_tcag/src/hmm; make; popd
pushd ${TOOLDIR}/ERDS-${ERDS}/erds_tcag/src; make; popd

cat <<EOF >${BINDIR}/erds_pipeline
#!/bin/bash
set -eu -o pipefail
export LC_ALL=en_US.UTF-8
perl ${TOOLDIR}/ERDS-${ERDS}/erds_tcag/src/erds_pipeline.pl \$@
EOF

chmod +x ${BINDIR}/erds_pipeline
