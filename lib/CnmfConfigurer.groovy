class CnmfConfigurer {
    static Map process(def params) {
        def cnmf = params?.cnmf
        def gpu = params.containsKey('cnmf_gpu') ? params.cnmf_gpu : [:]
        def container = ContainerConfigurer.process(params)
        boolean gpuEnabled = gpu?.enabled == true
        boolean gpuPrepareEnabled = gpuEnabled && gpu?.prepare?.enabled == true
        boolean gpuFactorizeEnabled = gpuEnabled && gpu?.factorize?.enabled == true

        Map runtime = [
            container: container,
            conda    : params?.rn_conda,
            scratch  : params?.rn_scratch
        ]

        [
            runtime   : runtime,
            common    : [
                k            : values(cnmf?.k),
                n_iter       : cnmf?.n_iter,
                n_workers    : cnmf?.n_workers,
                seed         : cnmf?.seed,
                local_density: cnmf?.local_density,
                save_h5ad    : cnmf?.save_h5ad
            ],
            preprocess: [
                n_variable: params?.preprocess?.n_variable
            ],
            prepare   : [
                engine   : gpuPrepareEnabled ? 'gpu' : 'cpu',
                label    : gpuPrepareEnabled ? (gpu?.prepare?.label ?: cnmf?.label) : cnmf?.label,
                scratch  : runtime.scratch,
                container: gpuPrepareEnabled
                    ? (gpu?.prepare?.container ?: gpu?.container ?: runtime.container)
                    : runtime.container,
                conda    : gpuPrepareEnabled
                    ? (gpu?.prepare?.conda ?: gpu?.conda ?: runtime.conda)
                    : runtime.conda,
                args     : gpuPrepareEnabled
                    ? joinArgs('--engine gpu', prepareArgs(gpu?.prepare))
                    : ''
            ],
            factorize : [
                engine          : gpuFactorizeEnabled ? 'gpu' : 'cpu',
                label           : gpuFactorizeEnabled
                    ? (gpu?.factorize?.label ?: gpu?.label ?: 'gpu_medium')
                    : cnmf?.label,
                scratch         : runtime.scratch,
                container       : gpuFactorizeEnabled
                    ? (gpu?.factorize?.container ?: gpu?.container ?: runtime.container)
                    : runtime.container,
                containerOptions: gpuFactorizeEnabled ? '--nv' : '',
                conda           : gpuFactorizeEnabled
                    ? (gpu?.factorize?.conda ?: gpu?.conda ?: runtime.conda)
                    : runtime.conda,
                args            : gpuFactorizeEnabled
                    ? joinArgs('--engine gpu', factorizeArgs(gpu, gpu?.factorize))
                    : ''
            ]
        ]
    }

    private static String prepareArgs(def prepare) {
        [
            option('--solver', prepare?.solver),
            option('--beta-loss', prepare?.beta_loss),
            option('--max-nmf-iter', prepare?.max_iter)
        ].findAll { it }.join(' ')
    }

    private static String factorizeArgs(def shared, def stageConfig) {
        def device = stageConfig?.device != null ? stageConfig.device : shared?.device
        def dtype = stageConfig?.dtype != null ? stageConfig.dtype : shared?.dtype
        def compileBlock = stageConfig?.compile_block != null ? stageConfig.compile_block : shared?.compile_block
        def allowTf32 = stageConfig?.allow_tf32 != null ? stageConfig.allow_tf32 : shared?.allow_tf32
        def compile = stageConfig?.compile != null ? stageConfig.compile : shared?.compile
        def eps = stageConfig?.eps != null ? stageConfig.eps : shared?.eps
        def checkEvery = stageConfig?.check_every != null ? stageConfig.check_every : shared?.check_every
        def batch = stageConfig?.batch != null ? stageConfig.batch : shared?.batch

        [
            option('--gpu-device', device),
            option('--gpu-dtype', dtype),
            option('--gpu-compile-block', compileBlock),
            booleanOption('--gpu-allow-tf32', allowTf32),
            booleanOption('--gpu-compile', compile),
            option('--gpu-eps', eps),
            option('--gpu-check-every', checkEvery),
            option('--gpu-batch', batch)
        ].findAll { it }.join(' ')
    }

    private static List<String> values(def value) {
        if (value == null) {
            return []
        }
        if (value instanceof Collection) {
            return value.collect { it.toString() }
        }
        value.toString().split(',').collect { it.trim() }.findAll { it }
    }

    private static String joinArgs(String... args) {
        args.findAll { it }.join(' ')
    }

    private static String option(String name, def value) {
        value != null ? "${name} ${value}" : ''
    }

    private static String booleanOption(String name, def value) {
        value != null && value.toString().trim().toLowerCase() in ['1', 'true', 'yes', 'on'] ? name : ''
    }
}
