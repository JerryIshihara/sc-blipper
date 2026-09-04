class ScviConfigurer {
    static Map process(def params) {
        def scvi = params?.preprocess?.scvi
        def sharedContainer = ContainerConfigurer.process(params)
        def container = scvi?.container ?: sharedContainer
        boolean gpuEnabled = scvi?.label?.toString()?.startsWith('gpu_')

        [
            label           : scvi?.label,
            container       : container,
            containerOptions: container && gpuEnabled ? '--nv' : '',
            conda           : scvi?.conda ?: params?.rn_conda
        ]
    }
}
