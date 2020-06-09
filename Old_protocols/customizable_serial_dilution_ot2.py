def get_values(*names):
    import json
    _all_values = json.loads("""{"pipette_type":"p300_single_gen2","dilution_factor":1.5,"num_of_dilutions":5,"total_mixing_volume":200,"tip_use_strategy":"never"}""")
    return [_all_values[n] for n in names]


metadata = {
    'protocolName': 'Customizable Serial Dilution',
    'author': 'Opentrons <protocols@opentrons.com>',
    'source': 'Protocol Library',
    'apiLevel': '2.2'
    }


def run(protocol_context):
    [pipette_type, dilution_factor, num_of_dilutions, total_mixing_volume,
        tip_use_strategy] = get_values( 
            'pipette_type', 'dilution_factor', 'num_of_dilutions',
            'total_mixing_volume', 'tip_use_strategy'
        )

    # labware
    trough = protocol_context.load_labware(
        'usascientific_12_reservoir_22ml', '2')
    liquid_trash = trough.wells()[0]
    

    plate = protocol_context.load_labware(
        'corning_96_wellplate_360ul_flat', '3')
    tipracks_200 = [
        protocol_context.load_labware('opentrons_96_filtertiprack_200ul', slot)
        for slot in ['1', '4']
    ]

    tipracks_20 = [
        protocol_context.load_labware('opentrons_96_filtertiprack_20ul', slot)
        for slot in ['5']
    ]

    tubes_rack = protocol_context.load_labware("opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap", 6)
    standard_tube = tubes_rack["A1"]

    pipette_p200 = protocol_context.load_instrument(
        pipette_type, mount='right', tip_racks=tipracks_200)

    pipette_p20 = protocol_context.load_instrument(
        "p20_single_gen2", mount='left', tip_racks=tipracks_20)


    transfer_volume = total_mixing_volume/dilution_factor
    diluent_volume = total_mixing_volume - transfer_volume
    standard_vol = 5

    
    ##Prefill plate with TE or whatever
    for col in plate.columns()[1:1+num_of_dilutions*2:2]:
        pipette_p200.distribute(
            diluent_volume, trough.wells()[0], [well for well in col][::2])

    ## Standards to plate 
    pipette_p200.distribute(
        diluent_volume-standard_vol, trough.wells()[0], [well for well in plate.columns()[0]][::2])
    pipette_p20.distribute(
        standard_vol, standard_tube, [well for well in plate.columns()[0]][::2])



    for row in plate.rows()[::2]:
        
        if tip_use_strategy == 'never':
            pipette_p200.pick_up_tip()


        for s, d in zip(row[:num_of_dilutions:2], row[1:1+num_of_dilutions*2:2]):

            pipette_p200.transfer(
                transfer_volume,
                s,
                d,
                mix_after=(3, total_mixing_volume/2),
                new_tip=tip_use_strategy
            )

            pipette_p200.transfer(
                transfer_volume,
                row[num_of_dilutions],
                liquid_trash,
                new_tip=tip_use_strategy,
                blow_out=True
            )


        if tip_use_strategy == 'never':
            pipette_p200.drop_tip()
