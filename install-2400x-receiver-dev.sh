#!/bin/sh
# Choose the receiver for an existing dump1090-fa install, then install the
# DEMOD_2400X build from this directory with install-2400x-dev.sh.
#
#   sudo sh install-2400x-receiver-dev.sh rtl              RTL-SDR dongle
#   sudo sh install-2400x-receiver-dev.sh rsp1b            SDRplay RSP1B, bias-T off
#   sudo sh install-2400x-receiver-dev.sh rsp1b-bias-t     SDRplay RSP1B, bias-T on
#                                                          (powers a masthead preamp)
#
# Only the receiver lines of /etc/default/dump1090-fa change (RECEIVER,
# RECEIVER_SERIAL, RECEIVER_GAIN, EXTRA_OPTIONS). The RTL-SDR settings are saved
# once as /etc/default/dump1090-fa.rtl and restored by 'rtl'.
#
# The RSP1B needs the SDRplay API (3.15) and the SoapySDRPlay3 module installed,
# and this directory built with SoapySDR support ("SoapySDR support: yes").
set -e

DEFAULTS=/etc/default/dump1090-fa
SAVED=$DEFAULTS.rtl
HERE="$(cd "$(dirname "$0")" && pwd)"

[ "$(id -u)" = 0 ] || { echo "run with sudo" >&2; exit 1; }
[ -f "$DEFAULTS" ] || { echo "$DEFAULTS missing, is the dump1090-fa package installed?" >&2; exit 1; }

# set KEY VALUE: replace the KEY= line in $DEFAULTS
set_key() {
    sed -i "s|^$1=.*|$1=$2|" "$DEFAULTS"
}
# saved KEY: the KEY= value from the saved RTL-SDR settings
saved() {
    sed -n "s|^$1=||p" "$SAVED"
}

case "$1" in
    rtl)
        if [ -f "$SAVED" ]; then
            for k in RECEIVER RECEIVER_SERIAL RECEIVER_GAIN EXTRA_OPTIONS; do
                set_key $k "$(saved $k)"
            done
        else
            set_key RECEIVER rtlsdr
            set_key RECEIVER_SERIAL ""
            set_key RECEIVER_GAIN ""
            set_key EXTRA_OPTIONS '""'
        fi
        ;;
    rsp1b|rsp1b-bias-t)
        "$HERE/dump1090" --help 2>&1 | grep -q -- --device-setting \
            || { echo "$HERE/dump1090 has no SoapySDR support, install libsoapysdr-dev and rebuild" >&2; exit 1; }
        if grep -q "^RECEIVER=rtlsdr" "$DEFAULTS"; then
            cp "$DEFAULTS" "$SAVED"
        fi
        set_key RECEIVER soapy
        set_key RECEIVER_SERIAL driver=sdrplay
        set_key RECEIVER_GAIN ""                  # maximum; adaptive gain works down from it
        if [ "$1" = rsp1b-bias-t ]; then
            set_key EXTRA_OPTIONS '"--device-setting biasT_ctrl=true"'
        else
            set_key EXTRA_OPTIONS '""'
        fi
        ;;
    *)
        echo "usage: sudo sh $0 rtl | rsp1b | rsp1b-bias-t" >&2
        exit 1
        ;;
esac

grep -E "^(RECEIVER|RECEIVER_SERIAL|RECEIVER_GAIN|EXTRA_OPTIONS)=" "$DEFAULTS"
sh "$HERE/install-2400x-dev.sh"
